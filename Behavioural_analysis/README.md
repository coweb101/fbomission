# Behavioural analysis

`FBOmission_Behavioural_Data_Analysis_anonym.R` analyses the behavioural data and the output of the PE modelling (folder `PE_modelling`). It contains five parts:

1. **Choice accuracy:** exclusion of participants performing below chance (binomial test) and a GLMM of choice accuracy (trial × learning context × reward probability), with plots.
2. **Posterior predictive check:** the same GLMM on choices simulated with the best-fitting reinforcement-learning model.
3. **Learning rates:** linear mixed model of the fitted learning rates of the best model, with plots.
4. **Parameter recovery:** correlation of fitted and recovered parameters, one scatter plot per parameter.
5. **Action values:** trial-wise action values of the correct and incorrect action per stimulus, with plots. This part also writes `aggregated_data/behav_and_pe_data_concatenated.csv`, which the EEG analysis (`EEG_multitemporal_analysis`) merges with the EEG data.

## Input data

All input files are in `aggregated_data/`. They contain only anonymous participant IDs (`subj_01` to `subj_48`).

| File | Content | Created by |
|---|---|---|
| `FBOmiss_behaviour_immediate_anonymized.csv` | Trial-wise behavioural data (48 participants × 480 trials) | `PE_modelling`, step 1 (`create_correct_response_key.m`) |
| `FBOmiss_learning_parameter.csv` | Fitted parameters of the best model (one row per participant) | `PE_modelling`, step 3 |
| `FBOmiss_immediate_Q_values_and_PEs.csv` | Trial-wise choice probabilities, action values and prediction errors of the best model | `PE_modelling`, step 3 |
| `FBOmiss_behav_recovered.csv` | Choices simulated with the best model (25 data sets per participant) | `PE_modelling`, step 4 |
| `FBOmiss_learning_parameter_recovered.csv` | Parameters recovered from the simulated data | `PE_modelling`, step 4 |

The last four files are copies of the files of the same name in `PE_modelling/comp_model_fit_export/`. If you rerun the modelling, copy them here again.

## How to run

1. Set R's working directory to this folder (`Behavioural_analysis`), e.g. with `setwd()` at the top of the script.
2. Create a folder `plots/` in this folder. The script saves all figures there.
3. Install the required packages (see below) and run the script from top to bottom.

`set.seed(17)` at the top of the script makes the jitter in the plots reproducible.

## Required R packages

data.table, plyr, dplyr, reshape2, tidyr, rstatix, lme4, lmerTest, buildmer, coin, emmeans, performance, psych, ggeffects, ggplot2, RColorBrewer, viridis

```r
install.packages(c("data.table", "plyr", "dplyr", "reshape2", "tidyr", "rstatix", "lme4",
                   "lmerTest", "buildmer", "coin", "emmeans", "performance", "psych",
                   "ggeffects", "ggplot2", "RColorBrewer", "viridis"))
```

## Output

- `plots/`: figures of choice accuracy (observed and simulated), learning rates, parameter recovery and action values.
- `aggregated_data/behav_and_pe_data_concatenated.csv`: behavioural data merged with the trial-wise action values and prediction errors, used by the EEG analysis.

## Note

The learning-rate analysis (part 3) and the parameter recovery plots (part 4) refer to the parameters of the model that was selected as the best model in `PE_modelling` (see `comp_model_fit_export/modelfits_immediate.csv` and the README in `PE_modelling`).
