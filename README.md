**FBOmission: Computational pipeline for quantifying valence and prediction error coding in omission-related responses measured via EEG in humans during reinforcement learning**

This repository stores code used to analyze behavioural and EEG data of a sample of 48 participants completing a probabilistic reinforcement learning task with omissions of feedback stimuli.
Importantly, these omissions entailed contextual valence based on the alternative possible outcome, which could either be a displayed reward or a displayed loss depending on the learning context (Get Reward vs Avoid Loss).

**Overview of pipeline**
1) PE_modelling: modelling of learning rates and single-trial prediction errors (PEs) based on behavioural choice data (Matlab)
2) EEG_preprocessing: preprocessing templates used for data preparation and cleaning (BrainVision Analyzer)
3) Behavioural_analysis: statistical analysis of choice accuracy data and modelled learning rates (R)
4) Multitemporal_analysis: statistical analysis of effects of contextual valence and PEs on display- and omission-related amplitudes (R)
