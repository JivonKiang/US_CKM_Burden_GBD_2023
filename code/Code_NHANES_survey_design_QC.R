# Complex-survey design smoke test for the NHANES CKM input layer.
# This script intentionally does not assign CKM stages. The candidate CVD
# indicator is an internal data-pipeline check and is not a Stage 4 estimate.

library(readr)
library(dplyr)
library(survey)

options(survey.lonely.psu = "adjust")

project <- Sys.getenv("CKM_PROJECT_ROOT", getwd())
input_path <- file.path(project, "data/NHANES/derived/nhanes_ckm_fasting_input_20260716.csv")
output_path <- file.path(project, "data/NHANES/derived/nhanes_survey_design_qc_20260716.csv")

data <- read_csv(input_path, show_col_types = FALSE) |>
  filter(cycle != "2021-2023 candidate", !is.na(pooled_2011_mar2020_fasting_weight)) |>
  mutate(
    survey_stratum = interaction(cycle, SDMVSTRA, drop = TRUE),
    survey_psu = interaction(cycle, SDMVPSU, drop = TRUE),
    # The exact Stage 4 definition remains eAppendix-dependent. This candidate
    # merely checks that the repeated heart failure/CHD/angina/MI/stroke fields
    # can be handled under the complex NHANES design.
    candidate_self_reported_cvd = as.integer(
      rowSums(across(c(MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F), ~ .x == 1), na.rm = TRUE) > 0
    )
  )

design <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_stratum,
  weights = ~pooled_2011_mar2020_fasting_weight,
  nest = TRUE,
  data = data
)

estimate <- svymean(~candidate_self_reported_cvd, design, na.rm = TRUE)
confidence <- confint(estimate)

result <- tibble(
  check = "candidate_self_reported_cvd_not_a_ckm_stage",
  analytic_n = nrow(data),
  design_degrees_of_freedom = degf(design),
  estimate_percent = 100 * coef(estimate)[[1]],
  ci_lower_percent = 100 * confidence[1, 1],
  ci_upper_percent = 100 * confidence[1, 2],
  weighting = "pooled fasting weights divided by 4.6; Taylor linearization",
  interpretation = "Internal survey-design check only. Do not report as CKM Stage 4."
)

write_csv(result, output_path)
print(result)
