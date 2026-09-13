# Local data inputs

This directory intentionally contains no third-party source exports or participant-level data.

## GBD inputs

Obtain the 24 annual GBD source CSV files directly from the appropriate IHME source and place them in `data/gbd/` (or set `CKM_GBD_DATA_DIR` to a different local directory). Use the filenames expected by `src/compute_gbd_component_trends.py` and verify that each file contains 34 annual rows from 1990 through 2023 with the columns `Year`, `Value`, `Lower bound`, and `Upper bound`.

The six component conditions overlap. Do not sum them or interpret them as CKM-stage prevalence. Follow current IHME terms and cite the original datasets as required.

## Optional NHANES input

The optional `src/nhanes_survey_design_qc.R` script requires a separately prepared, non-versioned local analytic input. Set `CKM_NHANES_QC_INPUT` to that local file and `CKM_NHANES_QC_OUTPUT` to a local output path. The script is an internal complex-survey design check only; it neither assigns CKM stages nor estimates CKM-stage prevalence.

No JAMA/NHANES table transcription, NHANES XPT file, participant-level input, derived analytic file, or WHO source export is distributed in this repository.
