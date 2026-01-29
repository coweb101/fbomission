### Pipeline for data preparation, analysis & visualization of prediction error
### and EEG data with FBOmission Data 
### (CW, 08/2024; 12/2025)

#### Routine & read data ####

remove(list = ls()) # clear workspace
getwd() # show current working directory
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss")

data <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_200800_anonym.csv") # wanna save data in between?

segment_length <- 250

aggr_data <- data
data_presentation <- subset(data, omission==-1)
data_omission <- subset(data, omission==1)

##### Multitemporal Analysis: Omission Data | Prediction Error | Frontocentral Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

# define which data set should be used
data <- data_omission
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

# code absolute PE variable
data$pe_absolute <- abs(data$PE)
# center around mean
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

# convert to data table (which is more efficient)
data <- data.table::as.data.table(data) 

# define electrodes that should be used in the analysis
electrodes <- c("F3", "Fz", "F4", "FC1", "FCz", "FC2")

# define number of iterations of the loop
segment_length <- 250

# set up empty data frame to save model results
pe_coefficients <- data.frame() 

## The following loop fits the same lmm with amplitudes of each of the 
## specified samplepoints and saves fixed effect coefficients, standard errors
## and p-values as well as follow-up tests for the interaction.

# (withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop))

set.seed(17) # for reproducibility

withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    downsampled_i <- seq(1, 1000, 4)[i] # current sample point
    # add current sample point to each electrode name 
    samples <- paste0(electrodes, "_", downsampled_i) 

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

# save data in working directory
data.table::fwrite(pe_coefficients,"aggregated_data/multitemp_pre200_post800_frontocentral_cluster_omission_trials.csv", row.names=F)


##### Multitemporal Analysis: Omission Data | Prediction Error | Centroparietal Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

# define which data set should be used
data <- data_omission
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

# code absolute PE variable
data$pe_absolute <- abs(data$PE)
# center around mean
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

# convert to data table (which is more efficient)
data <- data.table::as.data.table(data) 

# define electrodes that should be used in the analysis
electrodes <- c("CP1", "CP2", "P3", "Pz", "P4")

# define number of iterations of the loop
segment_length <- 250

# set up empty data frame to save model results
pe_coefficients <- data.frame() 

## The following loop fits the same lmm with amplitudes of each of the 
## specified samplepoints and saves fixed effect coefficients, standard errors
## and p-values as well as follow-up tests for the interaction.

# (withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop))

set.seed(17) # for reproducibility

withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    downsampled_i <- seq(1, 1000, 4)[i] # current sample point
    # add current sample point to each electrode name 
    samples <- paste0(electrodes, "_", downsampled_i) 
    
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

# save data in working directory
data.table::fwrite(pe_coefficients,"aggregated_data/multitemp_pre200_post800_centroparietal_cluster_omission_trials.csv", row.names=F)


##### Multitemporal Analysis: Display Data | Prediction Error | Frontocentral Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

# define which data set should be used
data <- data_presentation
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

# code absolute PE variable
data$pe_absolute <- abs(data$PE)
# center around mean
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

# convert to data table (which is more efficient)
data <- data.table::as.data.table(data) 

# define electrodes that should be used in the analysis
electrodes <- c("F3", "Fz", "F4", "FC1", "FCz", "FC2")

# define number of iterations of the loop
segment_length <- 250

# set up empty data frame to save model results
pe_coefficients <- data.frame() 

## The following loop fits the same lmm with amplitudes of each of the 
## specified samplepoints and saves fixed effect coefficients, standard errors
## and p-values as well as follow-up tests for the interaction.

# (withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop))

set.seed(17) # for reproducibility

withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    downsampled_i <- seq(1, 1000, 4)[i] # current sample point
    # add current sample point to each electrode name 
    samples <- paste0(electrodes, "_", downsampled_i) 
    
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

# save data in working directory
data.table::fwrite(pe_coefficients,"aggregated_data/multitemp_pre200_post800_frontocentral_cluster_display_trials.csv", row.names=F)


##### Multitemporal Analysis: Display Data | Prediction Error | Centroparietal Cluster #####

library(lmerTest) # generates p-values and automatically loads lme4

# define which data set should be used
data <- data_presentation
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

# code absolute PE variable
data$pe_absolute <- abs(data$PE)
# center around mean
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)

# convert to data table (which is more efficient)
data <- data.table::as.data.table(data) 

# define electrodes that should be used in the analysis
electrodes <- c("CP1", "CP2", "P3", "Pz", "P4")

# define number of iterations of the loop
segment_length <- 250

# set up empty data frame to save model results
pe_coefficients <- data.frame() 

## The following loop fits the same lmm with amplitudes of each of the 
## specified samplepoints and saves fixed effect coefficients, standard errors
## and p-values as well as follow-up tests for the interaction.

# (withCallingHandlers is wrapped around the for-loop to save non-convergence 
# warnings from lmer (which is not possible with warninigs() inside of the loop))

set.seed(17) # for reproducibility

withCallingHandlers({
  for (i in 1:segment_length) {
    
    ### Prepare data ###
    
    # define variables to subset
    id_variables <- c("id", "valence", "pe_absolute")
    
    downsampled_i <- seq(1, 1000, 4)[i] # current sample point
    # add current sample point to each electrode name 
    samples <- paste0(electrodes, "_", downsampled_i) 
    
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

# save data in working directory
data.table::fwrite(pe_coefficients,"aggregated_data/multitemp_pre200_post800_centroparietal_cluster_display_trials.csv", row.names=F)

