**FBOmission: Computational pipeline for quantifying valence and prediction error coding in omission-related responses measured via EEG in humans during reinforcement learning**

This repository stores code used to analyze behavioural and EEG data of a sample of 48 participants completing a probabilistic reinforcement learning task with omissions of feedback stimuli.
Importantly, these omissions entailed contextual valence based on the alternative possible outcome, which could either be a displayed reward or a displayed loss depending on the learning context (Get Reward vs Avoid Loss).

**Overview of pipeline**

The folders are meant to be run in this order; each folder has its own README with details.

1) `PE_modelling`: fitting of 42 reinforcement-learning models to the choice data, model comparison (BIC), simulation of trial-wise action values and prediction errors (PEs) with the two best models, parameter recovery and data simulation for posterior predictive checks (MATLAB)
2) `Behavioural_analysis`: statistical analysis of choice accuracy, posterior predictive checks, modelled learning rates, parameter recovery and action values; merges the behavioural data with the trial-wise PEs for the EEG analysis (R)
3) `EEG_multitemporal_analysis`: merging of EEG and PE data and sample-wise linear mixed models of the effects of contextual valence and PEs on display- and omission-related amplitudes, with plots (R)

**Data availability**

All data files in this repository contain only anonymous participant IDs (`subj_01` to `subj_48`). The raw logfiles are not public; the code that reads them is included for transparency. The preprocessed EEG data (read from .dat and .vmrk BVA export files, anonymized and merged with PE data of the two best models in `FBOMISS_Read_and_Merge_Data_700Pre_1500Post.R`; each of the two about 2.3 GB) are available at OSF via https://doi.org/10.17605/OSF.IO/FZCQJ.

**Models**

All analyses based on PEs are run for the two best models according to BIC. File names indicate the model:

| Suffix | Model (code) | Model (manuscript) | Learning rates |
|---|---|---|---|
| `_BIC` | 2b | M2.1 (Valence-Only model) | separate for positive and negative outcomes,<br>with update of the unchosen action |
| `_5b` | 5b | M5.1 (Valence-Context model) | separate for positive/negative × displayed/omitted outcomes,<br>with update of the unchosen action |

In modelfits.csv and the MATLAB code, the models are labelled by their features: M2.1 as “valence (+ unchosen update)”, M5.1 as “valence + appearance (+ unchosen update)”, where “appearance” refers to displayed vs omitted feedback.

**From scripts to figures and tables**

The scripts create the individual plots of all figures (except for the figure showing the task and topographies, see below); the final figures were arranged and annotated (e.g. panel labels) in PowerPoint. Figures 3–5 and S4 combine panels from the two PNG files created per multitemporal analysis (one with fixed effects and grand averages, one with simple slopes and grand averages by PE). Fig. S7 combines the grand averages by PE (bottom rows) of the first PNG file of each of the four analyses. All plots created from scripts are saved in a folder `plots/` within the respective analysis folder.

| Manuscript | Script | Output |
|:---|:---|:---|
| Fig. 1 (Task) | – (created in MS PowerPoint) | – |
| Fig. 2a, Table 1 (Choice Accuracy) | `Behavioural_analysis/`<br>`FBOmission_Behavioural_Data`<br>`_Analysis_anonym.R`, part 1 | `predicted_accuracy.png`, `behavioural_analysis_output.txt` |
| Fig. 2b/c, section 3.2 (Learning Rates) | `Behavioural_analysis/…`, part 3; <br> recovery of the learning-rate difference: part 4 | `learning_rates_predicted_BIC.png` /<br> `_5b.png`, `behavioural_analysis_output.txt` |
| Fig. 3 (Display, M2.1) | `EEG_multitemporal_analysis/FBOMISS_Multitemporal_Analysis`<br>`_Results_Plots_BIC.R` | `multitemp_results_display_…_BIC_<date>.png`,<br>`…_simple_slopes_BIC_<date>.png` |
| Fig. 4 (Omission, M2.1) | `…_Results_Plots_BIC.R` | `multitemp_results_omission_…_BIC_<date>.png`, `…_simple_slopes_BIC_<date>.png` |
| Fig. 5 (Omission, M5.1) | `…_Results_Plots_5b.R` | `multitemp_results_omission_…_5b_<date>.png`,<br>`…_simple_slopes_5b_<date>.png` |
| Tables 2 and 3 (Significant Time Windows), peak effects in sections 3.3.1 and 3.3.2 | `FBOMISS_Multitemporal_Analysis`<br>`_700Pre_1500Post.R` (LMMs), `FBOMISS_extraction_of_peak_effects.R` <br> (significant time windows and peak effects) | `multitemp results files/`, console output |
| Section 3.3.3 ((Omission) P3) | `FBOMISS_omission_P3_characterization.R` | console output |
| Table S1 (Model Comparison) | `PE_modelling/`<br>`FBOmiss_PE_modelling_pipeline.m` (step 2); `Behavioural_analysis/…`, part 6 | `modelfits.csv` |
| Fig. S1 (Parameter Recovery) | `Behavioural_analysis/…`, part 4 | `scatter_<parameter>_<model>.png` |
| Fig. S2, Table S2 (Posterior Predictive Check) | `Behavioural_analysis/…`, part 2 | `accuracy_posterior_prediction_BIC.png` /<br> `_5b.png`, `behavioural_analysis_output.txt` |
| Table S3 (Trial Counts) | `FBOMISS_Trial_Count_Table.R` | console output |
| Table S4 (Mapping) | – (created in MS Word) | – |
| Fig. S3, S5, S6 (Topographies) | – (created in BrainVision Analyzer) | – |
| Fig. S4 (Display, M5.1) | `…_Results_Plots_5b.R` | `multitemp_results_display_…_5b_<date>.png`,<br> `…_simple_slopes_5b_<date>.png` |
| Fig. S7 (ERPs by absolute PE) | `…_Results_Plots_BIC.R`, <br> `…_Results_Plots_5b.R` | bottom rows of<br>`multitemp_results_display_…_<model>_<date>.png`<br> and `multitemp_results_omission_…`<br>`_<model>_<date>.png` |
| Figure R1 (ramp subtraction, review only) | `EEG_multitemporal_analysis/`<br>`ramp_subtraction_rebuttal/FBOMISS_Follow_Up_Ramp_Subtraction.R` | `multitemp_results_omission_…_ramp_corrected_5b_<date>.png` |

**Citation**

If you use this code, please cite the publication:

Weber, C., Jocham, G., & Bellebaum, C. (under review). Cortical prediction error responses to reward and loss omission.

and the archived version of the code:

Weber, C., Jocham, G., & Bellebaum, C. (2026). *Computational pipeline for quantifying valence and prediction error coding in omission-related responses* (Version v[…]) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.[…]

**License**

The code in this repository is licensed under the MIT License (see `LICENSE`). The data on OSF are licensed under CC BY-NC 4.0 (Attribution-NonCommercial 4.0 International).
