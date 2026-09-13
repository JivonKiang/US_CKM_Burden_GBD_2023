# Optional complex-survey design smoke test for a local NHANES input.
# This script does not assign CKM stages or estimate a CKM-stage prevalence.

library(readr)
library(dplyr)
library(survey)

options(survey.lonely.psu = "adjust")

input_path <- Sys.getenv("CKM_NHANES_QC_INPUT")
output_path <- Sys.getenv("CKM_NHANES_QC_OUTPUT", unset = file.path("outputs", "nhanes_survey_design_qc.csv"))

if (input_path == "" || !file.exists(input_path)) {
  stop("Set CKM_NHANES_QC_INPUT to a local, non-versioned analytic input. No NHANES participant-level input is distributed here.")
}

data <- read_csv(input_path, show_col_types = FALSE) |>
  filter(cycle != "2021-2023 candidate", !is.na(pooled_2011_mar2020_fasting_weight)) |>
  mutate(
    survey_stratum = interaction(cycle, SDMVSTRA, drop = TRUE),
    survey_psu = interaction(cycle, SDMVPSU, drop = TRUE),
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
  weighting = "Caller-supplied pooled fasting weight; Taylor linearization",
  interpretation = "Internal survey-design check only. Do not report as CKM Stage 4."
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(result, output_path)
print(result)
