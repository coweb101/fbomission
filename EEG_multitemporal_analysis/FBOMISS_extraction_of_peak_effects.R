#### Extracting Peak Effects ####

### (CW, 06/2026)

##### Routine #####

remove(list = ls()) # clear workspace
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss") # set working directory

#### Define custom function for neighbour effects ####
fiveinarow_filter <- function(x) {
  x <- sort(unique(x))
  
  # Find break points in consecutive sequence
  breaks <- c(0, which(diff(x) != 1), length(x))
  
  # Split into runs of consecutive numbers
  runs <- lapply(seq_along(breaks[-1]),
                 function(i) x[(breaks[i] + 1):breaks[i + 1]])
  
  # Keep only runs of length >= 5
  out <- unlist(runs[lengths(runs) >= 5])
  return(out)
}

#### Display Frontocentral Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_display_trials_no_baseline.csv", quote="") # read data again

## Apply alpha correction
significanteffect_pe <- which(p.adjust(pe_coefficients$prob_pe, method = "BH") < .05)
significanteffect_valence <- which(p.adjust(pe_coefficients$prob_valence, method = "BH") < .05)
significanteffect_interaction <- which(p.adjust(pe_coefficients$prob_interaction, method = "BH") < .05)

## And only keep effects at 5 consecutive sample points

# apply on main effects and interaction
significanteffect_pe <- fiveinarow_filter(significanteffect_pe)
significanteffect_valence <- fiveinarow_filter(significanteffect_valence)
significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)

## For follow-up tests: alpha threshold adjustment based on number of significant 
## interactions multiplied with number of post-hoc tests (3)
significanteffect_pos_pe <- which(pe_coefficients$pvalue_pos  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_neg_pe <- which(pe_coefficients$pvalue_neg  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_diff_valence_pe <- which(pe_coefficients$pvalue_diff <  (.05/(length(significanteffect_interaction)*3)))

## Reduce post-hoc tests to tests where the interaction reached significance
significanteffect_pos_pe <- significanteffect_pos_pe[significanteffect_pos_pe%in%significanteffect_interaction]
significanteffect_neg_pe <- significanteffect_neg_pe[significanteffect_neg_pe%in%significanteffect_interaction]
significanteffect_diff_valence_pe <- significanteffect_diff_valence_pe[significanteffect_diff_valence_pe%in%significanteffect_interaction]

#### PE Main Effects ####

significanteffect_pe
diff(significanteffect_pe) # three time windows (one of them during first 200 ms)

# exclude those during choice review
significanteffect_pe <- significanteffect_pe[significanteffect_pe>50]

# view coefficients
pe_coefficients$coef_pe[significanteffect_pe]  # coefficients are all positively signed, therefore we look for the maximum

which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
# maximum at samplepoint 253
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))*4-704
# i.e. at 308 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[253, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[253, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)


#### FB Valence Main Effects ####

significanteffect_valence # effects broadly distributed
diff(significanteffect_valence)

## get more info for early time windows (before 400)
early_valence_effects <- significanteffect_valence[significanteffect_valence < 276]
pe_coefficients$coef_valence[early_valence_effects] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
# maximum at samplepoint 217
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))*4-704
#i.e. at 164 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[217, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[217, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

## get more info for later time windows

late_valence_effects <- significanteffect_valence[significanteffect_valence > 276]

pe_coefficients$coef_valence[late_valence_effects] # coefficients are first negatively signed only later positive again, therefore we look for the minimum at first

which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
# maximum at samplepoint 357
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))*4-704
# i.e. at 724 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[357, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[357, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


## get info for late positive effects
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
# maximum at samplepoint 257
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))*4-704
# i.e. at 724 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[548, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[548, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe


## simple slope of sign pos PEs
pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients after feedback onset are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# maximum at samplepoint 259
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# i.e. at 332 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[259, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[259, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

significanteffect_interaction_post <- significanteffect_interaction[significanteffect_interaction>13] # exclude those in the baseline

pe_coefficients$coef_neg[significanteffect_interaction_post] # coefficients after feedback onset are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))
# maximum at samplepoint 260
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))*4-704
# i.e. at 336 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[260, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[260, c("coef_neg", "se_neg", "pvalue_neg")],2)

## differencetest

pe_coefficients$coeff_diff[significanteffect_interaction_post] # coefficients after feedback onset are all negatively signed, therefore we look for the maximum
which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))
# maximum at samplepoint 262
which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))*4-704
# i.e. at 344 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[262, c("coeff_diff", "se_diff", "pvalue_diff")],3)
round(pe_coefficients[262, c("coeff_diff", "se_diff", "pvalue_diff")],2)


## feedback preceding period

## get simple slopes of PEs (or rather anticipatory effects before feedback onset)
pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients before feedback onset are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# maximum at samplepoint 1
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[1, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[1, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 2, i.e. at -696 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[2, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[2, c("coef_neg", "se_neg", "pvalue_neg")],2)

#### Display Centroparietal Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_display_trials_no_baseline.csv", quote="") # read data again

## Apply alpha correction
significanteffect_pe <- which(p.adjust(pe_coefficients$prob_pe, method = "BH") < .05)
significanteffect_valence <- which(p.adjust(pe_coefficients$prob_valence, method = "BH") < .05)
significanteffect_interaction <- which(p.adjust(pe_coefficients$prob_interaction, method = "BH") < .05)

## And only keep effects at 5 consecutive sample points

# apply on main effects and interaction
significanteffect_pe <- fiveinarow_filter(significanteffect_pe)
significanteffect_valence <- fiveinarow_filter(significanteffect_valence)
significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)

## For follow-up tests: alpha threshold adjustment based on number of significant 
## interactions multiplied with number of post-hoc tests (3)
significanteffect_pos_pe <- which(pe_coefficients$pvalue_pos  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_neg_pe <- which(pe_coefficients$pvalue_neg  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_diff_valence_pe <- which(pe_coefficients$pvalue_diff <  (.05/(length(significanteffect_interaction)*3)))

## Reduce post-hoc tests to tests where the interaction reached significance
significanteffect_pos_pe <- significanteffect_pos_pe[significanteffect_pos_pe%in%significanteffect_interaction]
significanteffect_neg_pe <- significanteffect_neg_pe[significanteffect_neg_pe%in%significanteffect_interaction]
significanteffect_diff_valence_pe <- significanteffect_diff_valence_pe[significanteffect_diff_valence_pe%in%significanteffect_interaction]

#### PE Main Effects ####

significanteffect_pe 
diff(significanteffect_pe) # five time windows, most of them again within the first 200ms
# exclude those 
significanteffect_pe <- significanteffect_pe[significanteffect_pe>55]

pe_coefficients$coef_pe[significanteffect_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
# maximum at samplepoint 265, i.e. at 356 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[265, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[265, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)


#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence)

pe_coefficients$coef_valence[significanteffect_valence] # vast majority of coefficients is positively signed, therefore we look at first for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
# maximum at samplepoint 254
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 312 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[254, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[254, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

# to get later, reversed effects, look for minimum
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))
# minimum at samplepoint 319
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 572 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[319, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[319, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe

## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# maximum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[1, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[1, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 5, i.e. at -684 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[5, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[5, c("coef_neg", "se_neg", "pvalue_neg")],2)




#### Omission Frontocentral Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_omission_trials_no_baseline.csv", quote="") # read data again

## Apply alpha correction
significanteffect_pe <- which(p.adjust(pe_coefficients$prob_pe, method = "BH") < .05)
significanteffect_valence <- which(p.adjust(pe_coefficients$prob_valence, method = "BH") < .05)
significanteffect_interaction <- which(p.adjust(pe_coefficients$prob_interaction, method = "BH") < .05)

## And only keep effects at 5 consecutive sample points



# apply on main effects and interaction
significanteffect_pe <- fiveinarow_filter(significanteffect_pe)
significanteffect_valence <- fiveinarow_filter(significanteffect_valence)
significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)

## For follow-up tests: alpha threshold adjustment based on number of significant 
## interactions multiplied with number of post-hoc tests (3)
significanteffect_pos_pe <- which(pe_coefficients$pvalue_pos  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_neg_pe <- which(pe_coefficients$pvalue_neg  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_diff_valence_pe <- which(pe_coefficients$pvalue_diff <  (.05/(length(significanteffect_interaction)*3)))

## Reduce post-hoc tests to tests where the interaction reached significance
significanteffect_pos_pe <- significanteffect_pos_pe[significanteffect_pos_pe%in%significanteffect_interaction]
significanteffect_neg_pe <- significanteffect_neg_pe[significanteffect_neg_pe%in%significanteffect_interaction]
significanteffect_diff_valence_pe <- significanteffect_diff_valence_pe[significanteffect_diff_valence_pe%in%significanteffect_interaction]

#### PE Main Effects ####

significanteffect_pe # not existent for the frontocentral cluster

#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence) # 4 distinct time windows (two early (before 400 ms) and two later ones (after 788 ms))

## get more info for early time windows

early_valence_effects <- significanteffect_valence[significanteffect_valence < 372]

pe_coefficients$coef_valence[early_valence_effects] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
# maximum at samplepoint 236
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))*4-704
#i.e. at 240 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[236, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[236, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

## get more info for later time windows

late_valence_effects <- significanteffect_valence[significanteffect_valence > 272]

pe_coefficients$coef_valence[late_valence_effects] # coefficients are all negatively signed, therefore we look for the minimum

which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
# maximum at samplepoint 374
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))*4-704
#i.e. at 792 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[374, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[374, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe


## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# maximum at samplepoint 2, i.e. at -696 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[2, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[2, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[1, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[1, c("coef_neg", "se_neg", "pvalue_neg")],2)

#### Omission Centroparietal Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_omission_trials_no_baseline.csv", quote="") # read data again

## Apply alpha correction
significanteffect_pe <- which(p.adjust(pe_coefficients$prob_pe, method = "BH") < .05)
significanteffect_valence <- which(p.adjust(pe_coefficients$prob_valence, method = "BH") < .05)
significanteffect_interaction <- which(p.adjust(pe_coefficients$prob_interaction, method = "BH") < .05)

## And only keep effects at 5 consecutive sample points

# apply on main effects and interaction
significanteffect_pe <- fiveinarow_filter(significanteffect_pe)
significanteffect_valence <- fiveinarow_filter(significanteffect_valence)
significanteffect_interaction <- fiveinarow_filter(significanteffect_interaction)

## For follow-up tests: alpha threshold adjustment based on number of significant 
## interactions multiplied with number of post-hoc tests (3)
significanteffect_pos_pe <- which(pe_coefficients$pvalue_pos  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_neg_pe <- which(pe_coefficients$pvalue_neg  < (.05/(length(significanteffect_interaction)*3)))
significanteffect_diff_valence_pe <- which(pe_coefficients$pvalue_diff <  (.05/(length(significanteffect_interaction)*3)))

## Reduce post-hoc tests to tests where the interaction reached significance
significanteffect_pos_pe <- significanteffect_pos_pe[significanteffect_pos_pe%in%significanteffect_interaction]
significanteffect_neg_pe <- significanteffect_neg_pe[significanteffect_neg_pe%in%significanteffect_interaction]
significanteffect_diff_valence_pe <- significanteffect_diff_valence_pe[significanteffect_diff_valence_pe%in%significanteffect_interaction]

#### PE Main Effects ####

significanteffect_pe 
diff(significanteffect_pe) # three time windows after 800 ms

pe_coefficients$coef_pe[significanteffect_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
# maximum at samplepoint 430
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))*4-704
# i.e. at 1016 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[430, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[430, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)


#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence) # one sustaining time window

pe_coefficients$coef_valence[significanteffect_valence] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
# maximum at samplepoint 238
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 252 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[238, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[238, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe

sum(significanteffect_pos_pe != significanteffect_neg_pe) # simple slopes at exact same time points significant

## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# maximum at samplepoint 2, i.e. at -696 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[2, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[2, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[1, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[1, c("coef_neg", "se_neg", "pvalue_neg")],2)





