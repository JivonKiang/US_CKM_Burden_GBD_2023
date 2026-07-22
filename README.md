# US Cardiovascular-Kidney-Metabolic Burden (GBD 2023)

Public code and aggregate source-data release for the manuscript:

> US Cardiovascular-Kidney-Metabolic Burden: Published NHANES Stages, Six Selected GBD Components, and WHO Risk Exposures

## Scope

This repository contains the non-identifiable, aggregate files needed to inspect the source-specific analyses. The GBD extracts cover the United States, both sexes, All Population, age-standardized rates, 1990-2023, and six non-additive components: ischemic heart disease, stroke, atrial fibrillation and flutter, lower-extremity peripheral arterial disease, diabetes mellitus, and chronic kidney disease. Four measures are supplied for each component: prevalence, incidence, deaths, and DALYs.

The repository does not contain individual-level NHANES records, NHANES XPT files, author contact information, submission forms, or the manuscript DOCX. The JAMA/NHANES files are published aggregate estimates transcribed from the cited article; the WHO file is an aggregate GHO extract.

## Data sources

Official IHME GHDx records and dataset DOIs:

- Cardiovascular burden estimates 1990-2023: https://ghdx.healthdata.org/record/ihme-data/cvd-1990-2023; DOI https://doi.org/10.6069/8HN3-3879
- CKD mortality, prevalence, and DALY estimates 1990-2023: https://ghdx.healthdata.org/record/ihme-data/gbd-2023-ckd-1990-2023; DOI https://doi.org/10.6069/BNCJ-5T64
- GBD 2023 risk exposure estimates 1990-2023: https://ghdx.healthdata.org/record/ihme-data/gbd-2023-risk-exposure-estimates-1990-2023; DOI https://doi.org/10.6069/SNHE-RS21

These are authoritative source-dataset identifiers. They are not the original GBD Compare query ID or receipt for the local CSV export. IHME data remain subject to the IHME free-of-charge non-commercial user agreement and terms of use. WHO data remain subject to the source system's terms.

## Layout

- `data/`: aggregate source data, the field-level GBD query manifest, and the derived trend table.
- `code/`: Python and R scripts used for data preparation, trend computation, figure construction, and survey-design quality control.
- `metadata/`: environment snapshots and official source-record notes.

## Reproducibility

The main aggregate analysis used Python 3.13.5 with pandas 2.3.3, NumPy 2.3.4, SciPy 1.16.3, statsmodels 0.14.5, and Matplotlib 3.10.7. Install the Python dependencies with `python -m pip install -r requirements.txt`. The scripts retain the source-specific analysis boundary: GBD and WHO estimates are not converted into NHANES CKM stages and are not mathematically collapsed into a single stage prevalence.

## Citation

Please cite the manuscript and the official IHME GHDx records listed above. A repository DOI will be added after this GitHub release is archived in Zenodo or another persistent repository.
