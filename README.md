**FBOmission: Computational pipeline for quantifying valence and prediction error coding in omission-related responses measured via EEG in humans during reinforcement learning**

This repository stores code used to analyze behavioural and EEG data of a sample of 48 participants completing a probabilistic reinforcement learning task with omissions of feedback stimuli.
Importantly, these omissions entailed contextual valence based on the alternative possible outcome, which could either be a displayed reward or a displayed loss depending on the learning context (Get Reward vs Avoid Loss).

**Overview of pipeline**

The folders are meant to be run in this order; each folder has its own README with details.

1) `PE_modelling`: fitting of 42 reinforcement-learning models to the choice data, model comparison (BIC), simulation of trial-wise action values and prediction errors (PEs) with the two best models, parameter recovery and posterior predictive checks (MATLAB)
2) `Behavioural_analysis`: statistical analysis of choice accuracy, posterior predictive checks, modelled learning rates, parameter recovery and action values; merges the behavioural data with the trial-wise PEs for the EEG analysis (R)
3) `EEG_multitemporal_analysis`: merging of EEG and PE data and sample-wise linear mixed models of the effects of contextual valence and PEs on display- and omission-related amplitudes, with plots (R)

**Models**

All analyses based on PEs are run for the two best models according to BIC. File names indicate the model:

| Suffix | Model (code) | Model (manuscript) | Learning rates |
|---|---|---|---|
| `_BIC` | 2b | M2.1 | separate for positive and negative outcomes, with update of the unchosen action |
| `_5b` | 5b | M5.1 | separate for positive/negative × presented/omitted outcomes, with update of the unchosen action |

**Data availability**

All data files in this repository contain only anonymous participant IDs (`subj_01` to `subj_48`). The raw logfiles are not public; the code that reads them is included for transparency. The preprocessed EEG data (read from .dat and .vmrk BVA export files, anonymized and merged with PE data in `FBOMISS_Read_and_Merge_Data_700Pre_1500Post.R`; about 1.7 GB) are available at https://osf.io/fzcqj/.
