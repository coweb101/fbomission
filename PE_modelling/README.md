# PE modelling

`FBOmiss_PE_modelling_pipeline.m` runs the complete PE modelling pipeline. It calls the scripts and functions in the subfolders, so the folder structure and names must not be changed, and MATLAB's current folder must be `PE_modelling`.

The pipeline has four steps:

1. Reading the logfiles and preparing the behavioural data
2. Fitting 36 reinforcement-learning models and comparing them (−LL, AIC, BIC)
3. Simulating trial-wise action values and prediction errors with the best model by BIC (model 2b)
4. Parameter recovery and posterior predictive check with simulated data

Step 1 needs the raw logfiles, which are not public; its code is included for transparency. Steps 2–4 can be reproduced with the anonymized data in `interim_datasets`. To do so, run `rng(17, 'twister')` and then the pipeline from section `%% 2.` onwards (step 1 uses no random numbers, so the results are identical to a full run). Steps 1 and 2 take about 8 hours with 50 fitting iterations per model.

Requirements: MATLAB with the Optimization Toolbox (`fmincon`) and the Statistics and Machine Learning Toolbox.

**Note on file paths:** The MATLAB scripts were written and run under Windows and use backslashes in file paths (e.g. `save('comp_model_fit_export\fit_model_1_ll', ...)`). On macOS or Linux, replace these backslashes with forward slashes (`/`).
