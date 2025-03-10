# Description
This repository contains data documentation and code for the analysis in the manuscript titled "Feasibility of applying pharmacoepidemiologic drug-drug interaction screening methods to nursing home residents: An application to clopidogrel."
## Repository Contents
- `data_documentation/` - Contains files describing the data sources, key variables, and broad cohort creation steps, and technical notes on methods.>.
- `code/` - The programs used for data analysis.
- `LICENSE` - The license under which this repository is shared.
- `README.md` - This file, providing an overview of the repository.
## Data Documentation
The `data_documentation/` directory contains the following files:
- `1_Metadata_nhddi_screening_github.xlsx` - contains the list of input datasets and years of data used in the analysis; broad cohort creation steps, description of key variables in the final and near-final analytic datasets; list of drugs within drug classess of interest (described more in technical notes).
- `2_Technical_notes_nhddi_screening_github.xlsx` - includes additional technical notes about study methods that were not included in the published manuscript or supplaments. Also includes additional figures that help illustrate some of the cohort creation rules and processes. 
## Code
The `code/` directory contains the following programs:
- `1_conditional_poisson_models_nhddi_screening_github.sas` - includes essential code used to run conditional poisson models. Assumes person-day level analytic dataset has already been constructed. 
- `2_semibayes_adjustment_nhddi_screening_github.sas` - includes essential code to apply semi-Bayes adjustment to reduce false postivies from multiple estimation. Uses estimate and variance outputs from program 1. Adjustment is applied separately for set of estimates for each outcome.

Programs were run in sequence to produce the study findings.
Cohort creation programs have not been included, though a broad description of these steps can be found in the data documentation section. 

Additional information (and code) for identifying drug-observable time within the nursing home population can be found in the upcoming publication from Harris et al. "Identifying observable medication time for US nursing home residents using Medicare claims: A tutorial and case study"

