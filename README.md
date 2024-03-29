# Quantifying imperfect recall in the POEM score

This repository contains the code developed for the article by [Duverdier et al. (2024), "Evaluation of measurement errors in the Patient-Oriented Eczema Measure (POEM) outcome"](https://doi.org/10.1111/cea.14441), published in Clinical & Experimental Allergy.

The Patient-Oriented Eczema Measure (POEM) is the recommended the core outcome measure of eczema symptoms perceived by patients in clinical trials and practice. POEM is reported by recalling the presence/absence of symptoms that occurred in the last seven days. The FDA has highlighted the importance of considering the recall period of patient-reported outcome measures. This project investigated measurement errors in the Patient-Oriented Eczema Measure (POEM) score due to the imperfect recall of symptoms.

## Requirements/Installation

The code is written in the R language for statistical computing (version 4.2.0) and the probabilistic programming language [Stan](https://mc-stan.org/) for the models. 
Package dependencies can be found in [`renv.lock`](renv.lock).

## File structure
Files for the analysis conducted in this project are:
 - [`01_missing_values_imputation.R`](01_missing_values_imputation.R) : Impute missing daily symptom presence/absence using a [Markov chain model](https://ghurault.github.io/EczemaPred/articles/MC.html) implemented in the [EczemaPred](https://ghurault.github.io/EczemaPred/index.html) R package.
 - [`02_recall_error.R`](02_recall_error.R) : Compare d-POEM (derived from daily diaries) and r-POEM (recalled) scores, and calculate recall bias and recall noise in the POEM score.
 - [`03a_recall_model_check.R`](03a_recall_model_check.R) : Prior predictive check for the recalled days model.
 - [`03b_recall_model_fit.R`](03b_recall_model_fit.R) : Fit the recalled days model on the data and plot the results of the model.
 
Common functions used throughout the project can be found in [`functions.R`](functions.R).

The code for the Stan recalled days model developed in this project can be found in [`Model/recalled_days_model.stan`](Model/recalled_days_model.stan).

## License
The open source version of this project is licensed under the GPLv3 license, which can be seen in the [LICENSE](LICENSE) file.
