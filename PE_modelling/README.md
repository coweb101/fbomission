The file FBOmiss_PE_modelling_pipeline.m was run to read logfiles, prepare raw data and fit and compare 
the different models to the behavioural data. It works as a batch file that calls functions in the subfolders 
within this folder. The folder structure (and names) should therefore not be changed.

The first step runs on the non-public raw data files and was published for transparency.
Everything from step 2 onwards can be reproduced using the anonymized files provided in the folder interim_datasets.

**Note on file paths:**
MATLAB scripts were written and run under Windows and use backslashes in file paths (e.g. `save('comp_model_fit_export\fit_model_1_ll', ...)`).
On macOS or Linux, replace these backslashes with forward slashes (`/`).
