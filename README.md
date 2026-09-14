# US Cardiovascular-Kidney-Metabolic Burden (GBD 2023)

Reproducibility code and source documentation for a source-respecting descriptive analysis of cardiovascular-kidney-metabolic (CKM) burden and monitoring in the United States.

## Scope

This is a code-only release. It contains no third-party source exports, participant-level data, NHANES XPT files, derived analytic inputs, manuscript files, author-contribution records, credentials, or transcribed JAMA/NHANES tables. The six GBD disease components can overlap and must not be summed or interpreted as a CKM-stage prevalence measure. The optional NHANES script is a survey-design smoke test for non-staging checks only; it does not assign CKM stages.

## Data access

The GBD input CSV files and their permissions remain the responsibility of the user. Obtain data directly from the appropriate source, review its current terms, and place the 24 annual source files in `data/gbd/` using the names expected by `src/compute_gbd_component_trends.py`. The script expects 34 annual observations from 1990 through 2023 and columns `Year`, `Value`, `Lower bound`, and `Upper bound`.

IHME data downloadable from IHME websites may be used, shared, modified, or built upon by non-commercial users under the IHME Free-of-Charge Non-Commercial User Agreement. The data are not redistributed or relicensed by this repository. Reuse of JAMA Network tables, figures, or selected text requires permission through RightsLink, so no transcribed JAMA/NHANES table is distributed here. Users who independently obtain WHO data must follow the applicable dataset terms, retain prescribed attribution, and avoid implying WHO endorsement.

The optional NHANES survey-design check requires a separately prepared local analytic input. It is intentionally excluded from this repository. Set `CKM_NHANES_QC_INPUT` and `CKM_NHANES_QC_OUTPUT` before running it.

## Run

```text
python -m pip install -r requirements.txt
python src/compute_gbd_component_trends.py
```

Set `CKM_GBD_DATA_DIR` to an alternative GBD input directory if needed. The command writes `outputs/gbd_component_trends_1990_2023.csv` by default; set `CKM_OUTPUT_DIR` to change the output directory.

For the optional NHANES design check:

```text
Rscript src/nhanes_survey_design_qc.R
```

## Methods boundary

The primary trend estimand is the average annual percentage change from a BIC-selected continuous two-segment log-linear model. Candidate breakpoints are prespecified interior years; the selected model uses one-lag HAC covariance. Intervals do not include breakpoint-selection uncertainty or annual GBD input uncertainty. The global log-linear result is a sensitivity diagnostic.

This repository provides code for descriptive trend summaries and source documentation. It is not a CKM-stage estimator, a cross-source validation framework, or a causal analysis.

## Citation

Please cite the associated manuscript when it is published. Cite source datasets according to the requirements of their original providers. The relevant GBD records are documented in `metadata/official_source_records.md`.

## License

The `LICENSE` file applies only to original code. Third-party source-data terms remain controlling.
