# EEG multitemporal analysis

Sample-wise linear mixed models (LMMs) of the effects of contextual valence and the absolute prediction error (|PE|) on feedback-locked EEG amplitudes, separately for displayed and omitted feedback, and the corresponding plots. All PE-based analyses are run for the two best models from `PE_modelling` (`_BIC` = M2.1, model 2b; `_5b` = M5.1, model 5b).

## Scripts (in the order in which they are run)

| Script | What it does | Input | Output |
|---|---|---|---|
| `FBOMISS_Read_and_Merge_Data_700Pre_1500Post.R` | Reads the single-trial EEG exports (BrainVision Analyzer), adds empty rows for trials rejected during preprocessing, checks the trial order against the behavioural data via the feedback markers and merges EEG and PE data | BrainVision exports (.dat and .vmrk) and ID key (created within the MATLAB pipeline; both not public), `behav_and_pe_data_concatenated_BIC.csv` / `_5b.csv` (from `Behavioural_analysis`) | `eeg_and_pe_data_merged_700pre_1500post_anonym_BIC.csv` / `_5b.csv` |
| `FBOMISS_Trial_Count_Table.R` | Number of valid trials per condition | merged data | console output |
| `FBOMISS_Multitemporal_Analysis_700Pre_1500Post.R` | One LMM per sample point, separately for display and omission trials and for the frontocentral and centroparietal cluster, for both models | merged data (`_BIC`, `_5b`) | `multitemp_pre700_post1500_<cluster>_cluster_<condition>_trials_no_baseline_<model>.csv` |
| `FBOMISS_Multitemporal_Analysis_Results_Plots_BIC.R`, `..._5b.R` | Plots of the LMM coefficients (with significant effects after Benjamini-Hochberg correction and at least five consecutive sample points) and grand averages (with 40 ms moving average) | LMM results and merged data of the respective model | figures in `plots/` |
| `FBOMISS_extraction_of_peak_effects.R` | Extracts the peak effects of the significant time windows | LMM results | console output |
| `FBOMISS_omission_P3_characterization.R` | Characterizes the P3-like waveform in omission (and display) trials | merged data (_BIC) | console output (latencies and peaks) |

## Analysis details

- Segments from −700 to 1500 ms relative to feedback onset, 250 Hz (550 sample points; sample point *i* corresponds to *i* × 4 − 704 ms), no baseline correction.
- Electrode clusters: frontocentral (F3, Fz, F4, FC1, FCz, FC2) and centroparietal (CP1, CP2, P3, Pz, P4).
- Model per sample point: `amplitude ~ |PE| × valence + (1 + |PE| + |PE|:valence | participant) + (1 | electrode)`, with |PE| rescaled to [−1, 1] and centred and valence centred; simple slopes of |PE| for positive and negative outcomes with emmeans.
- Participants whose |PE| hardly varies (SD of |PE| < .05 in display or omission trials) are excluded from all PE analyses, as their participant-specific PE slopes cannot be estimated. The criterion is applied to the PEs of each model separately.
- Grand averages for the PE plots use terciles of |PE| formed within participant and outcome valence.

## Result files

`multitemp results files/` contains the LMM results (one row per sample point; estimates, standard errors, degrees of freedom, t and p values for each effect and the simple slopes). The scripts read and write data in `aggregated_data/` relative to the working directory; to rerun the plots, copy the result files there.

## Data

The merged EEG and PE data (`eeg_and_pe_data_merged_700pre_1500post_anonym_BIC.csv`, `…_5b.csv`) are too large for GitHub and are available at [OSF](https://doi.org/10.17605/OSF.IO/FZCQJ). To rerun the analyses, place them in `aggregated_data/`. `FBOMISS_Read_and_Merge_Data_700Pre_1500Post.R` created these files from the single-trial BrainVision exports and is included for transparency only. It needs the non-anonymized exports and the participant key, which are not public.

## How to run

Set the working directory at the top of each script (`setwd()`) to the folder that contains `aggregated_data/` and create a folder `plots/`.

## Required R packages

data.table, plyr, lmerTest (loads lme4), emmeans, effectsize, scales, viridis, yarrr

```r
install.packages(c("data.table", "plyr", "lmerTest", "emmeans", "effectsize",
                   "scales", "viridis", "yarrr"))
```
