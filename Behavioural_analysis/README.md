# Behavioural analysis

`FBOmission_Behavioural_Data_Analysis_anonym.R` analyses the behavioural data and the output of the PE modelling (folder `PE_modelling`). All model-based parts are run for the two best models according to BIC: files ending in `_BIC` refer to M2.1 (model 2b, lowest BIC), files ending in `_5b` to M5.1 (model 5b, second-lowest BIC).

The script contains six parts:

1. **Choice accuracy:** exclusion of participants performing below chance (binomial test) and a GLMM of choice accuracy (trial × learning context × reward probability), with plots.
2. **Posterior predictive check:** the same GLMM on choices simulated with each of the two models (25 simulated data sets per participant), with plots.
3. **Learning rates:** linear mixed models of the fitted (cube-root transformed) learning rates of both models (M2.1: valence; M5.1: valence × feedback appearance), with plots.
4. **Parameter recovery:** correlation of fitted and recovered parameters of both models (recovered parameters averaged across the 25 simulated data sets), one scatter plot per parameter.
5. **Action values:** trial-wise action values of the correct and incorrect action per stimulus, with plots. This part also writes the behavioural data merged with the trial-wise action values and prediction errors, which the EEG analysis (`EEG_multitemporal_analysis`) merges with the EEG data.
6. **Model fits:** table of the fit indices of all 42 models (for Supplementary Table S1).

## Input data

All input files are in `aggregated_data/`. They contain only anonymous participant IDs (`subj_01` to `subj_48`).

| File | Content | Created by |
|---|---|---|
| `FBOmiss_behaviour_immediate_anonymized.csv` | Trial-wise behavioural data (48 participants × 480 trials) | `PE_modelling`, step 1 (`create_correct_response_key.m`) |
| `modelfits.csv` | Fit indices (−LL, BIC, AIC) of all 42 models | `PE_modelling`, step 2 |
| `FBOmiss_learning_parameter_BIC.csv`, `_5b.csv` | Fitted parameters (one row per participant) | `PE_modelling`, step 3 |
| `FBOmiss_immediate_Q_values_and_PEs_BIC.csv`, `_5b.csv` | Trial-wise choice probabilities, action values and prediction errors | `PE_modelling`, step 3 |
| `FBOmiss_behav_recovered_BIC.csv`, `_5b.csv` | Choices simulated with the model (25 data sets per participant) | `PE_modelling`, step 4 |
| `FBOmiss_learning_parameter_recovered_BIC.csv`, `_5b.csv` | Parameters recovered from the simulated data | `PE_modelling`, step 4 |

All files except the first are copies of the files of the same name in `PE_modelling/comp_model_fit_export/`. If you rerun the modelling, copy them here again.

## How to run

1. Set R's working directory to this folder (`Behavioural_analysis`), e.g. with `setwd()` at the top of the script.
2. Create a folder `plots/` in this folder. The script saves all figures there.
3. Install the required packages (see below) and run the script from top to bottom.

`set.seed(17)` at the top of the script makes the jitter in the plots reproducible.

## Required R packages

data.table, plyr, reshape2, lme4, lmerTest, emmeans, performance, psych, ggeffects, viridis, knitr

```r
install.packages(c("data.table", "plyr", "reshape2", "lme4", "lmerTest", "emmeans",
                   "performance", "psych", "ggeffects", "viridis", "knitr"))
```

## Output

- `plots/`: choice accuracy (`predicted_accuracy.png`), posterior predictive checks (`accuracy_posterior_prediction_<model>.png`), learning rates (`learning_rates_predicted_<model>.png`), parameter recovery (`scatter_<parameter>_<model>.png`) and action values (`action_values_correct_<model>.png`, `action_values_full_range_<model>.png`), with `<model>` = `BIC` or `5b`.
- `aggregated_data/behav_and_pe_data_concatenated_BIC.csv` and `aggregated_data/behav_and_pe_data_concatenated_5b.csv`: behavioural data merged with the trial-wise action values and prediction errors of M2.1 and M5.1, used by the EEG analysis.

## Note

The model comparison and the selection of the two models are described in the README of `PE_modelling` (see `comp_model_fit_export/modelfits.csv`).
