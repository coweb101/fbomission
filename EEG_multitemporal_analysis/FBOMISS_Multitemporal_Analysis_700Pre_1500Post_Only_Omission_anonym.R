#### FBOmission: Multitemporal Analysis | Extended Segment (-700 to 1500 ms) ####

### (CW, 02/2025)

##### Routine & read data #####

remove(list = ls()) # clear workspace
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss") # set working directory
data <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_7001500_anonym.csv", quote="") # read concatenated data

data_omission <- subset(data, omission==1)
data_presentation <- subset(data, omission==-1)
data_all <- data

segment_length <- 550 
set.seed(17) # for reproducibility


##### Multitemporal Analysis: Omission Data | Prediction Error | Frontocentral electrodes #####

library(lmerTest) # generates p-values and automatically loads lme4

# subset data
#aggr_data <- data # save full data in other object
data <- data_omission # subset only omission trials
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

gc()

# Take absolute values of PEs (as we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

data <- data.table::as.data.table(data) # convert to data table which seems to be more efficient

# define variable electrodes
electrodes <- c("F3", "Fz", "F4", "FC1", "FCz", "FC2")

pe_coefficients <- data.frame() # set up empty data frame to save coefficients

## The following loop fits the same lmm with amplitudes of each of the exported 
## samplepoints and saves coefficients (involving the effect of prediction 
## error), standard errors and p-values.

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)
# info: runs ~33min with example data

# caution: in the current form data must be of type data.table because of indexing style in loop
withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    # add current sample point to each electrode name
    # (data is already downsampled)
    samples <- paste0(electrodes, "_", i) 
    
    # concatenate variable names
    all_variables <- c(id_variables, samples)
    
    # subset these variables
    data_temp <- data[,..all_variables] # the dots are needed here because I 
    # converted the data frame to a data table which is allegedly more time 
    # efficient
    
    # melt data so that all amplitudes are in one column
    data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                  variable.name = "electrode",
                                  value.name = "samplepointx")
    
    # ensure that amplitude data is numeric
    data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
    
    
    
    ### Calculate mixed model ###
    
    # DV: amplitude of sample point i of all electrodes specified before loop
    # Fixed Effects: absolute prediction error & feedback valence (within-subjects)
    # and their interaction
    # Random effects: intercept for each participant and electrode, 
    # by-participant slopes for PE and the interaction between PE & Valence
    
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*valence + 
                                               (1+ pe_absolute + pe_absolute:valence|id) + (1|electrode),
                                             REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    
    # save singular fit warning
    if (isSingular(samplepoint_regression)) {
      pe_coefficients[i, "singular"] <- 1
    } else {pe_coefficients[i, "singular"] <- 0}
    
    pe_coefficients[i, "convergence_false"] <- 0
    # enter default value for non-convergence (which is overwritten when a 
    # warning is caught with calling handlers (see below)
    
    
    
    ### Extract fixed effects coefficients, SEs, dfs, t- and p-values ###
    
    # Intercept
    pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "prob_intercept")] <- 
      coef(summary(samplepoint_regression))["(Intercept)",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # PE
    pe_coefficients[i,c("coef_pe", "se_pe", "df_pe", "t_pe", "prob_pe")] <- 
      coef(summary(samplepoint_regression))["pe_absolute",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # FB Valence
    pe_coefficients[i,c("coef_valence", "se_valence", "df_valence", "t_valence", "prob_valence")] <- 
      coef(summary(samplepoint_regression))["valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # Interaction PE x FB Valence
    pe_coefficients[i,c("coef_interaction_pe_valence", "se_interaction_pe_valence",
                        "df_interaction_pe_valence", "t_interaction_pe_valence",
                        "prob_interaction_pe_valence")] <- 
      coef(summary(samplepoint_regression))["pe_absolute:valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    
    
    ### Conduct simple slope tests for positive and negative PEs ###
    emmeans <- emmeans::emtrends(samplepoint_regression,
                                 specs=pairwise~valence+pe_absolute, 
                                 var="pe_absolute",
                                 infer=TRUE,
                                 adj="none",
                                 lmer.df="asymp")
    
    emtrends <- as.data.frame(summary(emmeans)$emtrends) # extract simple slope tests
    diff <- as.data.frame(summary(emmeans)$contrasts) # extract test of difference of slopes
    
    # Extract marginal effect coefficients, SEs, and p-values for positive and
    # negative PEs
    pe_coefficients[i,"coef_neg"] <- emtrends[which(emtrends$valence==-1),c("pe_absolute.trend")]
    pe_coefficients[i,"coef_pos"] <- emtrends[which(emtrends$valence==1),c("pe_absolute.trend")]
    
    pe_coefficients[i,"se_neg"] <- emtrends[which(emtrends$valence==-1),c("SE")]
    pe_coefficients[i,"se_pos"] <- emtrends[which(emtrends$valence==1),c("SE")]
    
    pe_coefficients[i,"pvalue_neg"] <- emtrends[which(emtrends$valence==-1),c("p.value")]
    pe_coefficients[i,"pvalue_pos"] <- emtrends[which(emtrends$valence==1),c("p.value")]
    
    # & of difference test between the simple slopes
    pe_coefficients[i,"coeff_diff"] <- diff[1,c("estimate")]
    pe_coefficients[i,"se_diff"] <- diff[1,c("SE")]
    pe_coefficients[i,"pvalue_diff"] <- diff[1,c("p.value")]
    
    ### Clean up & monitor ###
    
    # remove temporary objects to free memory
    rm("emmeans", "emtrends", "data_temp") 
    print(i) # print number of current sample point to monitor progress
    
  }
}, warning = function(w){
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})


# add number of samplepoints for plotting
pe_coefficients$time <- 1:nrow(pe_coefficients) 

#### Save data in working directory
pe_coefficients_electrode_cluster <- pe_coefficients
data.table::fwrite(pe_coefficients_electrode_cluster, "aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_omission_trials.csv", row.names=F)



##### Multitemporal Analysis: Omission Data | Prediction Error | Centroparietal electrodes #####

library(lmerTest) # generates p-values and automatically loads lme4

# subset data
#aggr_data <- data # save full data in other object
data <- data_omission # subset only omission trials
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

gc()

# Take absolute values of PEs (as we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

data <- data.table::as.data.table(data) # convert to data table which seems to be more efficient

# define variable electrodes
electrodes <- c("CP1", "CP2", "P3", "Pz", "P4")

pe_coefficients <- data.frame() # set up empty data frame to save coefficients

## The following loop fits the same lmm with amplitudes of each of the exported 
## samplepoints and saves coefficients (involving the effect of prediction 
## error), standard errors and p-values.

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)
# info: runs ~33min with example data

# caution: in the current form data must be of type data.table because of indexing style in loop
withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    # add current sample point to each electrode name
    # (data is already downsampled)
    samples <- paste0(electrodes, "_", i) 
    
    # concatenate variable names
    all_variables <- c(id_variables, samples)
    
    # subset these variables
    data_temp <- data[,..all_variables] # the dots are needed here because I 
    # converted the data frame to a data table which is allegedly more time 
    # efficient
    
    # melt data so that all amplitudes are in one column
    data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                  variable.name = "electrode",
                                  value.name = "samplepointx")
    
    # ensure that amplitude data is numeric
    data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
    
    
    
    ### Calculate mixed model ###
    
    # DV: amplitude of sample point i of all electrodes specified before loop
    # Fixed Effects: absolute prediction error & feedback valence (within-subjects)
    # and their interaction
    # Random effects: intercept for each participant and electrode, 
    # by-participant slopes for PE and the interaction between PE & Valence
    
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*valence + 
                                               (1+ pe_absolute + pe_absolute:valence|id) + (1|electrode),
                                             REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    
    # save singular fit warning
    if (isSingular(samplepoint_regression)) {
      pe_coefficients[i, "singular"] <- 1
    } else {pe_coefficients[i, "singular"] <- 0}
    
    pe_coefficients[i, "convergence_false"] <- 0
    # enter default value for non-convergence (which is overwritten when a 
    # warning is caught with calling handlers (see below)
    
    
    
    ### Extract fixed effects coefficients, SEs, dfs, t- and p-values ###
    
    # Intercept
    pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "prob_intercept")] <- 
      coef(summary(samplepoint_regression))["(Intercept)",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # PE
    pe_coefficients[i,c("coef_pe", "se_pe", "df_pe", "t_pe", "prob_pe")] <- 
      coef(summary(samplepoint_regression))["pe_absolute",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # FB Valence
    pe_coefficients[i,c("coef_valence", "se_valence", "df_valence", "t_valence", "prob_valence")] <- 
      coef(summary(samplepoint_regression))["valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # Interaction PE x FB Valence
    pe_coefficients[i,c("coef_interaction_pe_valence", "se_interaction_pe_valence",
                        "df_interaction_pe_valence", "t_interaction_pe_valence",
                        "prob_interaction_pe_valence")] <- 
      coef(summary(samplepoint_regression))["pe_absolute:valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    
    
    ### Conduct simple slope tests for positive and negative PEs ###
    emmeans <- emmeans::emtrends(samplepoint_regression,
                                 specs=pairwise~valence+pe_absolute, 
                                 var="pe_absolute",
                                 infer=TRUE,
                                 adj="none",
                                 lmer.df="asymp")
    
    emtrends <- as.data.frame(summary(emmeans)$emtrends) # extract simple slope tests
    diff <- as.data.frame(summary(emmeans)$contrasts) # extract test of difference of slopes
    
    # Extract marginal effect coefficients, SEs, and p-values for positive and
    # negative PEs
    pe_coefficients[i,"coef_neg"] <- emtrends[which(emtrends$valence==-1),c("pe_absolute.trend")]
    pe_coefficients[i,"coef_pos"] <- emtrends[which(emtrends$valence==1),c("pe_absolute.trend")]
    
    pe_coefficients[i,"se_neg"] <- emtrends[which(emtrends$valence==-1),c("SE")]
    pe_coefficients[i,"se_pos"] <- emtrends[which(emtrends$valence==1),c("SE")]
    
    pe_coefficients[i,"pvalue_neg"] <- emtrends[which(emtrends$valence==-1),c("p.value")]
    pe_coefficients[i,"pvalue_pos"] <- emtrends[which(emtrends$valence==1),c("p.value")]
    
    # & of difference test between the simple slopes
    pe_coefficients[i,"coeff_diff"] <- diff[1,c("estimate")]
    pe_coefficients[i,"se_diff"] <- diff[1,c("SE")]
    pe_coefficients[i,"pvalue_diff"] <- diff[1,c("p.value")]
    
    ### Clean up & monitor ###
    
    # remove temporary objects to free memory
    rm("emmeans", "emtrends", "data_temp") 
    print(i) # print number of current sample point to monitor progress
    
  }
}, warning = function(w){
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})


# add number of samplepoints for plotting
pe_coefficients$time <- 1:nrow(pe_coefficients) 

#### Save data in working directory
pe_coefficients_electrode_cluster <- pe_coefficients
data.table::fwrite(pe_coefficients_electrode_cluster, "aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_omission_trials.csv", row.names=F)





##### Multitemporal Analysis: Display Data | Prediction Error | Frontocentral electrodes #####

library(lmerTest) # generates p-values and automatically loads lme4

# subset data
#aggr_data <- data # save full data in other object
data <- data_presentation # subset only omission trials
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

gc()

# Take absolute values of PEs (as we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

data <- data.table::as.data.table(data) # convert to data table which seems to be more efficient

# define variable electrodes
electrodes <- c("F3", "Fz", "F4", "FC1", "FCz", "FC2")

pe_coefficients <- data.frame() # set up empty data frame to save coefficients

## The following loop fits the same lmm with amplitudes of each of the exported 
## samplepoints and saves coefficients (involving the effect of prediction 
## error), standard errors and p-values.

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)
# info: runs ~33min with example data

# caution: in the current form data must be of type data.table because of indexing style in loop

withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    # add current sample point to each electrode name
    # (data is already downsampled)
    samples <- paste0(electrodes, "_", i) 
    
    # concatenate variable names
    all_variables <- c(id_variables, samples)
    
    # subset these variables
    data_temp <- data[,..all_variables] # the dots are needed here because I 
    # converted the data frame to a data table which is allegedly more time 
    # efficient
    
    # melt data so that all amplitudes are in one column
    data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                  variable.name = "electrode",
                                  value.name = "samplepointx")
    
    # ensure that amplitude data is numeric
    data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
    
    
    
    ### Calculate mixed model ###
    
    # DV: amplitude of sample point i of all electrodes specified before loop
    # Fixed Effects: absolute prediction error & feedback valence (within-subjects)
    # and their interaction
    # Random effects: intercept for each participant and electrode, 
    # by-participant slopes for PE and the interaction between PE & Valence
    
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*valence + 
                                               (1+ pe_absolute + pe_absolute:valence|id) + (1|electrode),
                                             REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    
    # save singular fit warning
    if (isSingular(samplepoint_regression)) {
      pe_coefficients[i, "singular"] <- 1
    } else {pe_coefficients[i, "singular"] <- 0}
    
    pe_coefficients[i, "convergence_false"] <- 0
    # enter default value for non-convergence (which is overwritten when a 
    # warning is caught with calling handlers (see below)
    
    
    
    ### Extract fixed effects coefficients, SEs, dfs, t- and p-values ###
    
    # Intercept
    pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "prob_intercept")] <- 
      coef(summary(samplepoint_regression))["(Intercept)",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # PE
    pe_coefficients[i,c("coef_pe", "se_pe", "df_pe", "t_pe", "prob_pe")] <- 
      coef(summary(samplepoint_regression))["pe_absolute",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # FB Valence
    pe_coefficients[i,c("coef_valence", "se_valence", "df_valence", "t_valence", "prob_valence")] <- 
      coef(summary(samplepoint_regression))["valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # Interaction PE x FB Valence
    pe_coefficients[i,c("coef_interaction_pe_valence", "se_interaction_pe_valence",
                        "df_interaction_pe_valence", "t_interaction_pe_valence",
                        "prob_interaction_pe_valence")] <- 
      coef(summary(samplepoint_regression))["pe_absolute:valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    
    
    ### Conduct simple slope tests for positive and negative PEs ###
    emmeans <- emmeans::emtrends(samplepoint_regression,
                                 specs=pairwise~valence+pe_absolute, 
                                 var="pe_absolute",
                                 infer=TRUE,
                                 adj="none",
                                 lmer.df="asymp")
    
    emtrends <- as.data.frame(summary(emmeans)$emtrends) # extract simple slope tests
    diff <- as.data.frame(summary(emmeans)$contrasts) # extract test of difference of slopes
    
    # Extract marginal effect coefficients, SEs, and p-values for positive and
    # negative PEs
    pe_coefficients[i,"coef_neg"] <- emtrends[which(emtrends$valence==-1),c("pe_absolute.trend")]
    pe_coefficients[i,"coef_pos"] <- emtrends[which(emtrends$valence==1),c("pe_absolute.trend")]
    
    pe_coefficients[i,"se_neg"] <- emtrends[which(emtrends$valence==-1),c("SE")]
    pe_coefficients[i,"se_pos"] <- emtrends[which(emtrends$valence==1),c("SE")]
    
    pe_coefficients[i,"pvalue_neg"] <- emtrends[which(emtrends$valence==-1),c("p.value")]
    pe_coefficients[i,"pvalue_pos"] <- emtrends[which(emtrends$valence==1),c("p.value")]
    
    # & of difference test between the simple slopes
    pe_coefficients[i,"coeff_diff"] <- diff[1,c("estimate")]
    pe_coefficients[i,"se_diff"] <- diff[1,c("SE")]
    pe_coefficients[i,"pvalue_diff"] <- diff[1,c("p.value")]
    
    ### Clean up & monitor ###
    
    # remove temporary objects to free memory
    rm("emmeans", "emtrends", "data_temp") 
    print(i) # print number of current sample point to monitor progress
    
  }
}, warning = function(w){
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})



## Calculate conditional coefficients of the effect of prediction error (i.e.
## separate coefficients for each feedback timing group and valence (which are
## dummy-coded and therefore additive))

#pe_coefficients$coef_pos <- pe_coefficients$coef_pe
####pe_coefficients$coef_pos_immediate <- pe_coefficients$coef_pe + pe_coefficients$coef_pe_absolute_fb_timing_dummy
#pe_coefficients$coef_neg <- pe_coefficients$coef_pe + pe_coefficients$coef_pe_absolute_valence_dummy
####pe_coefficients$coef_neg_immediate <- pe_coefficients$coef_pe + pe_coefficients$coef_pe_absolute_fb_timing_dummy + 
####  pe_coefficients$coef_pe_absolute_valence_dummy + pe_coefficients$coef_pe_absolute_fb_timing_dummy_valence_dummy

#pe_coefficients$time <- 1:segment_length # add number of samplepoints for plotting
# as my model loop above was interrupted due to 
pe_coefficients$time <- 1:nrow(pe_coefficients) 

#### Save data in working directory
pe_coefficients_electrode_cluster <- pe_coefficients
data.table::fwrite(pe_coefficients_electrode_cluster, "aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_display_trials.csv", row.names=F)



##### Multitemporal Analysis: Display Data | Prediction Error | Centroparietal electrodes #####

library(lmerTest) # generates p-values and automatically loads lme4

# subset data
#aggr_data <- data # save full data in other object
data <- data_presentation # subset only omission trials
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

gc()

# Take absolute values of PEs (as we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

data <- data.table::as.data.table(data) # convert to data table which seems to be more efficient

# define variable electrodes
electrodes <- c("CP1", "CP2", "P3", "Pz", "P4")

pe_coefficients <- data.frame() # set up empty data frame to save coefficients

## The following loop fits the same lmm with amplitudes of each of the exported 
## samplepoints and saves coefficients (involving the effect of prediction 
## error), standard errors and p-values.

# withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop)
# info: runs ~33min with example data

# caution: in the current form data must be of type data.table because of indexing style in loop

withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    # add current sample point to each electrode name
    # (data is already downsampled)
    samples <- paste0(electrodes, "_", i) 
    
    # concatenate variable names
    all_variables <- c(id_variables, samples)
    
    # subset these variables
    data_temp <- data[,..all_variables] # the dots are needed here because I 
    # converted the data frame to a data table which is allegedly more time 
    # efficient
    
    # melt data so that all amplitudes are in one column
    data_temp <- data.table::melt(data_temp, id.vars=id_variables,
                                  variable.name = "electrode",
                                  value.name = "samplepointx")
    
    # ensure that amplitude data is numeric
    data_temp$samplepointx <- as.numeric(data_temp$samplepointx)
    
    
    
    ### Calculate mixed model ###
    
    # DV: amplitude of sample point i of all electrodes specified before loop
    # Fixed Effects: absolute prediction error & feedback valence (within-subjects)
    # and their interaction
    # Random effects: intercept for each participant and electrode, 
    # by-participant slopes for PE and the interaction between PE & Valence
    
    samplepoint_regression <- lmerTest::lmer(samplepointx ~ pe_absolute*valence + 
                                               (1+ pe_absolute + pe_absolute:valence|id) + (1|electrode),
                                             REML=T, data = data_temp, control = lmerControl(optimizer="bobyqa"))
    
    # save singular fit warning
    if (isSingular(samplepoint_regression)) {
      pe_coefficients[i, "singular"] <- 1
    } else {pe_coefficients[i, "singular"] <- 0}
    
    pe_coefficients[i, "convergence_false"] <- 0
    # enter default value for non-convergence (which is overwritten when a 
    # warning is caught with calling handlers (see below)
    
    
    
    ### Extract fixed effects coefficients, SEs, dfs, t- and p-values ###
    
    # Intercept
    pe_coefficients[i,c("coef_intercept", "se_intercept", "df_intercept", "t_intercept", "prob_intercept")] <- 
      coef(summary(samplepoint_regression))["(Intercept)",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # PE
    pe_coefficients[i,c("coef_pe", "se_pe", "df_pe", "t_pe", "prob_pe")] <- 
      coef(summary(samplepoint_regression))["pe_absolute",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # FB Valence
    pe_coefficients[i,c("coef_valence", "se_valence", "df_valence", "t_valence", "prob_valence")] <- 
      coef(summary(samplepoint_regression))["valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    # Interaction PE x FB Valence
    pe_coefficients[i,c("coef_interaction_pe_valence", "se_interaction_pe_valence",
                        "df_interaction_pe_valence", "t_interaction_pe_valence",
                        "prob_interaction_pe_valence")] <- 
      coef(summary(samplepoint_regression))["pe_absolute:valence",
                                            c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
    
    
    
    ### Conduct simple slope tests for positive and negative PEs ###
    emmeans <- emmeans::emtrends(samplepoint_regression,
                                 specs=pairwise~valence+pe_absolute, 
                                 var="pe_absolute",
                                 infer=TRUE,
                                 adj="none",
                                 lmer.df="asymp")
    
    emtrends <- as.data.frame(summary(emmeans)$emtrends) # extract simple slope tests
    diff <- as.data.frame(summary(emmeans)$contrasts) # extract test of difference of slopes
    
    # Extract marginal effect coefficients, SEs, and p-values for positive and
    # negative PEs
    pe_coefficients[i,"coef_neg"] <- emtrends[which(emtrends$valence==-1),c("pe_absolute.trend")]
    pe_coefficients[i,"coef_pos"] <- emtrends[which(emtrends$valence==1),c("pe_absolute.trend")]
    
    pe_coefficients[i,"se_neg"] <- emtrends[which(emtrends$valence==-1),c("SE")]
    pe_coefficients[i,"se_pos"] <- emtrends[which(emtrends$valence==1),c("SE")]
    
    pe_coefficients[i,"pvalue_neg"] <- emtrends[which(emtrends$valence==-1),c("p.value")]
    pe_coefficients[i,"pvalue_pos"] <- emtrends[which(emtrends$valence==1),c("p.value")]
    
    # & of difference test between the simple slopes
    pe_coefficients[i,"coeff_diff"] <- diff[1,c("estimate")]
    pe_coefficients[i,"se_diff"] <- diff[1,c("SE")]
    pe_coefficients[i,"pvalue_diff"] <- diff[1,c("p.value")]
    
    ### Clean up & monitor ###
    
    # remove temporary objects to free memory
    rm("emmeans", "emtrends", "data_temp") 
    print(i) # print number of current sample point to monitor progress
    
  }
}, warning = function(w){
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message # save convergence message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})


# add number of samplepoints for plotting
pe_coefficients$time <- 1:nrow(pe_coefficients) 

#### Save data in working directory
pe_coefficients_electrode_cluster <- pe_coefficients
data.table::fwrite(pe_coefficients_electrode_cluster, "aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_display_trials.csv", row.names=F)





