#### Extracting Peak Effects Model 2b (best BIC) ####

### (CW, 10/2026)

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

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_display_trials_no_baseline_BIC.csv", quote="") # read data again

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


#### Time Windows for Table for Displayed Feedback - FC cluster (M2.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe
diff(significanteffect_pe) # three time windows (one of them during first 200 ms)

split(significanteffect_pe, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))


# exclude those during choice review
significanteffect_pe <- significanteffect_pe[significanteffect_pe>50]

# view coefficients
pe_coefficients$coef_pe[significanteffect_pe]  # coefficients are all positively signed, therefore we look for the maximum

which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
# maximum at samplepoint 257
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))*4-704
# i.e. at 324 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)


#### FB Valence Main Effects ####

significanteffect_valence # effects broadly distributed
diff(significanteffect_valence)

split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))


## get more info for early time windows (before 400)
early_valence_effects <- significanteffect_valence[significanteffect_valence < 276]
pe_coefficients$coef_valence[early_valence_effects] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
# maximum at samplepoint 217
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))*4-704
#i.e. at 164 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

## get more info for later time windows

late_valence_effects <- significanteffect_valence[significanteffect_valence > 275]

pe_coefficients$coef_valence[late_valence_effects] # coefficients are first negatively signed only later positive again, therefore we look for the minimum at first

which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
# minimum at samplepoint 357
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))*4-704
# i.e. at 724 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


## get info for late positive effects
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
# maximum at samplepoint 548
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))*4-704
# i.e. at 1488 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe

split(significanteffect_interaction, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))

split(significanteffect_pos_pe, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))

split(significanteffect_neg_pe, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))

split(significanteffect_diff_valence_pe, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


## simple slope of sign pos PEs
pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients after feedback onset are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# maximum at samplepoint 258
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# i.e. at 328 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

significanteffect_interaction_post <- significanteffect_interaction[significanteffect_interaction>50] # exclude those in the choice review period

pe_coefficients$coef_neg[significanteffect_interaction_post] # coefficients after feedback onset are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))
# maximum at samplepoint 278
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))*4-704
# i.e. at 408 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))
, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))
, c("coef_neg", "se_neg", "pvalue_neg")],2)

## differencetest

pe_coefficients$coeff_diff[significanteffect_interaction_post] # coefficients after feedback onset are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))
# minimum at samplepoint 258
which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))*4-704
# i.e. at 328 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))
, c("coeff_diff", "se_diff", "pvalue_diff")],3)
round(pe_coefficients[which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))
, c("coeff_diff", "se_diff", "pvalue_diff")],2)


## feedback preceding period

## get simple slopes of PEs (or rather anticipatory effects before feedback onset)
pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients before feedback onset are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# minimum at samplepoint 1
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 4, i.e. at -688 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
, c("coef_neg", "se_neg", "pvalue_neg")],2)

#### Display Centroparietal Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_display_trials_no_baseline_BIC.csv", quote="") # read data again

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

#### Time Windows for Table for Displayed Feedback - CP cluster (M2.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe 
diff(significanteffect_pe) # five time windows, many of them again within the first 200ms or in the feedback preceding period

split(significanteffect_pe, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))


# exclude those 
significanteffect_pe <- significanteffect_pe[significanteffect_pe>126]

pe_coefficients$coef_pe[significanteffect_pe] # first few very small coefficients negatively signed, others are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))*4-704
# maximum at samplepoint 267, i.e. at 364 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
, c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)


#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence)
split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))


pe_coefficients$coef_valence[significanteffect_valence] # vast majority of coefficients is positively signed, therefore we look at first for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
# maximum at samplepoint 254
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 312 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

# to get later, reversed effects, look for minimum
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))
# minimum at samplepoint 418
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 968 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe

split(significanteffect_interaction, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))

split(significanteffect_pos_pe, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))

split(significanteffect_neg_pe, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))

split(significanteffect_diff_valence_pe, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))



## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are mostly negatively signed, therefore we look for the minimum at first
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# minimum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],2)

# look for later max
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# maximum at samplepoint 264, i.e. at 352 ms relative to fb omission onset (samplepoint*4-704)

round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],2)


## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
, c("coef_neg", "se_neg", "pvalue_neg")],2)




#### Omission Frontocentral Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_omission_trials_no_baseline_BIC.csv", quote="") # read data again

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

#### Time Windows for Table for Omitted Feedback - FC cluster (M2.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe # not existent for the frontocentral cluster

#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence) 
split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
# 6 distinct time windows, some before feedback omission onset, most afterwards

### get more info for early time windows

early_valence_effects <- significanteffect_valence[significanteffect_valence < 176]

pe_coefficients$coef_valence[early_valence_effects] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
# maximum at samplepoint 135
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))*4-704
#i.e. at -164 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

## get more info for later time windows

late_valence_effects <- significanteffect_valence[significanteffect_valence > 175]

pe_coefficients$coef_valence[late_valence_effects] # coefficients are mostly positvely, late ones negatively signed, therefore we look for the maximum at first

which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
# maximum at samplepoint 269
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))*4-704
#i.e. at 372 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

# then look for the minimum

which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
# min at samplepoint 375
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))*4-704
#i.e. at 796 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

#### Interactions ####

significanteffect_interaction
diff(significanteffect_interaction) 
split(significanteffect_interaction, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
# one significant time window between -700 and -592

significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))



## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# minimum at samplepoint 2, i.e. at -696 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe])), c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe])), c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 2, i.e. at -696 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe])), c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe])), c("coef_neg", "se_neg", "pvalue_neg")],2)

#### Omission Centroparietal Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_omission_trials_no_baseline_BIC.csv", quote="") # read data again

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

#### Time Windows for Table for Omitted Feedback - CP cluster (M2.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe # no effects

#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence)
split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
# 4 time windows, all well after omission onset

pe_coefficients$coef_valence[significanteffect_valence] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
# maximum at samplepoint 237
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 244 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
, c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe

## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# minimum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
, c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
, c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
, c("coef_neg", "se_neg", "pvalue_neg")],2)






#### Extracting Peak Effects Model 5b (second best BIC) ####

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

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_display_trials_no_baseline_5b.csv", quote="") # read data again

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

#### Time Windows for Table for Displayed Feedback - FC cluster (M5.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe
diff(significanteffect_pe) 
split(significanteffect_pe, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
# one time window between 284 and 360 ms

# view coefficients
pe_coefficients$coef_pe[significanteffect_pe]  # coefficients are all positively signed, therefore we look for the maximum

which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
# maximum at samplepoint 260 (2b: 257)
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))*4-704
# i.e. at 336 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)


#### FB Valence Main Effects ####

significanteffect_valence # effects broadly distributed
diff(significanteffect_valence)

split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))


## get more info for early time windows (before 400)
early_valence_effects <- significanteffect_valence[significanteffect_valence < 276]
pe_coefficients$coef_valence[early_valence_effects] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
# maximum at samplepoint 217 (as for 2b)
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))*4-704
#i.e. at 164 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

## get more info for later time windows

late_valence_effects <- significanteffect_valence[significanteffect_valence > 276]

pe_coefficients$coef_valence[late_valence_effects] # coefficients are first negatively signed only later positive again, therefore we look for the minimum at first

which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
# minimum at samplepoint 358 (2b: 357)
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))*4-704
# i.e. at 728 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


## get info for late positive effects
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
# maximum at samplepoint 548 (as for 2b)
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))*4-704
# i.e. at 1488 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe


## simple slope of sign pos PEs
pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients after feedback onset are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# maximum at samplepoint 258 (as for 2b)
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# i.e. at 328 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

significanteffect_interaction_post <- significanteffect_interaction[significanteffect_interaction>50] # exclude those in the choice review period

pe_coefficients$coef_neg[significanteffect_interaction_post] # coefficients partly positive, partly negative, therefore we look at first for the maximum

which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))
# maximum at samplepoint 253 (for model2b: 278)
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))*4-704
# i.e. at 308 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_interaction_post]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],2)

# look for min
which(pe_coefficients$coef_neg == min(pe_coefficients$coef_neg[significanteffect_interaction_post]))
# min at samplepoint 269
which(pe_coefficients$coef_neg == min(pe_coefficients$coef_neg[significanteffect_interaction_post]))*4-704
# i.e. at 372 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == min(pe_coefficients$coef_neg[significanteffect_interaction_post]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == min(pe_coefficients$coef_neg[significanteffect_interaction_post]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],2)



## differencetest

pe_coefficients$coeff_diff[significanteffect_interaction_post] # coefficients after feedback onset are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))
# minimum at samplepoint 258 (as for 2b)
which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))*4-704
# i.e. at 328 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))
                      , c("coeff_diff", "se_diff", "pvalue_diff")],3)
round(pe_coefficients[which(pe_coefficients$coeff_diff == min(pe_coefficients$coeff_diff[significanteffect_interaction_post]))
                      , c("coeff_diff", "se_diff", "pvalue_diff")],2)


## feedback preceding period

## get simple slopes of PEs (or rather anticipatory effects before feedback onset)
pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients before feedback onset are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# min at samplepoint 1 (as for 2b)
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # no effects for neg PEs

#### Display Centroparietal Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_display_trials_no_baseline_5b.csv", quote="") # read data again

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

#### Time Windows for Table for Displayed Feedback - CP cluster (M5.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe # no main effects

#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence)
split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))


pe_coefficients$coef_valence[significanteffect_valence] # vast majority of coefficients is positively signed, therefore we look at first for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
# maximum at samplepoint 254 (as for 2b)
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 312 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

# to get later, reversed effects, look for minimum
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))
# minimum at samplepoint 319 (2b: 418)
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 572 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[significanteffect_valence]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)


#### Interactions ####

significanteffect_interaction
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe

## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are mostly negatively signed, therefore we look for the minimum at first
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# minimum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],2)

# look for later max
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# maximum at samplepoint 265 (2b: 264), i.e. at 352 ms relative to fb omission onset (samplepoint*4-704)

round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == max(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],2)


## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 5, i.e. at -684 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],2)




#### Omission Frontocentral Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_omission_trials_no_baseline_5b.csv", quote="") # read data again

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

#### Time Windows for Table for Omitted Feedback - FC cluster (M5.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe # not existent for the frontocentral cluster

#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence) 
split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
# 6 distinct time windows, some before feedback omission onset, most afterwards

### get more info for early time windows

early_valence_effects <- significanteffect_valence[significanteffect_valence < 176]

pe_coefficients$coef_valence[early_valence_effects] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
# maximum at samplepoint 135 (as for 2b)
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))*4-704
#i.e. at -164 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[early_valence_effects]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

## get more info for later time windows

late_valence_effects <- significanteffect_valence[significanteffect_valence > 176]

pe_coefficients$coef_valence[late_valence_effects] # coefficients are mostly positvely, late ones negatively signed, therefore we look for the maximum at first

which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))
# maximum at samplepoint 269 (as for 2b)
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects]))*4-704
#i.e. at 372 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

# then look for the minimum

which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))
# min at samplepoint 375 (as for 2b)
which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects]))*4-704
#i.e. at 796 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == min(pe_coefficients$coef_valence[late_valence_effects])), c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

#### Interactions ####

significanteffect_interaction
diff(significanteffect_interaction) 
split(significanteffect_interaction, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
# one significant time window between -700 and -592 (as for 2b)

significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))



## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
# minimum at samplepoint 16, i.e. at -640 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe])), c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe])), c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe])), c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe])), c("coef_neg", "se_neg", "pvalue_neg")],2)

#### Omission Centroparietal Cluster ####

pe_coefficients <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_omission_trials_no_baseline_5b.csv", quote="") # read data again

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

#### Time Windows for Table for Omitted Feedback - CP cluster (M5.1) ####
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_interaction*4-704, cumsum(c(TRUE, diff(significanteffect_interaction) != 1)))
split(significanteffect_pos_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pos_pe) != 1)))
split(significanteffect_neg_pe*4-704, cumsum(c(TRUE, diff(significanteffect_neg_pe) != 1)))
split(significanteffect_diff_valence_pe*4-704, cumsum(c(TRUE, diff(significanteffect_diff_valence_pe) != 1)))


#### PE Main Effects ####

significanteffect_pe 
split(significanteffect_pe, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
split(significanteffect_pe*4-704, cumsum(c(TRUE, diff(significanteffect_pe) != 1)))
# 7 time windows, first starts shortly before omission onset

pe_coefficients$coef_pe[significanteffect_pe] # coefficients are at first negatively, then positively signed, therefore we look for the maximum at first

which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))
# maximum at samplepoint 381
which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe]))*4-704
# i.e. at 820 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[which(pe_coefficients$coef_pe == max(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)

# then look for minimum

which(pe_coefficients$coef_pe == min(pe_coefficients$coef_pe[significanteffect_pe]))
# minimum at samplepoint 181
which(pe_coefficients$coef_pe == min(pe_coefficients$coef_pe[significanteffect_pe]))*4-704
# i.e. at 20 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pe == min(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],3)
round(pe_coefficients[which(pe_coefficients$coef_pe == min(pe_coefficients$coef_pe[significanteffect_pe])), c("coef_pe", "se_pe", "df_pe", "t_pe" , "prob_pe", "effectsize_pe")],2)


#### FB Valence Main Effects ####

significanteffect_valence 

diff(significanteffect_valence)
split(significanteffect_valence, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
split(significanteffect_valence*4-704, cumsum(c(TRUE, diff(significanteffect_valence) != 1)))
# 3 time windows

pe_coefficients$coef_valence[significanteffect_valence] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
# maximum at samplepoint 238 (2b: 237)
which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))*4-704
#i.e. at 248 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],3)
round(pe_coefficients[which(pe_coefficients$coef_valence == max(pe_coefficients$coef_valence[significanteffect_valence]))
                      , c("coef_valence", "se_valence", "df_valence", "t_valence" , "prob_valence", "effectsize_valence")],2)

#### Interactions ####

significanteffect_interaction
significanteffect_interaction*4-704
significanteffect_pos_pe
significanteffect_neg_pe
significanteffect_diff_valence_pe

## simple slope of pos PEs

pe_coefficients$coef_pos[significanteffect_pos_pe] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))*4-704
# minimum at samplepoint 5, i.e. at -684 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],2)

## simple slope of neg PEs

pe_coefficients$coef_neg[significanteffect_neg_pe] # coefficients are all positively signed, therefore we look for the maximum
which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
# maximum at samplepoint 1, i.e. at -700 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],3)
round(pe_coefficients[which(pe_coefficients$coef_neg == max(pe_coefficients$coef_neg[significanteffect_neg_pe]))
                      , c("coef_neg", "se_neg", "pvalue_neg")],2)


# later interaction effect


significanteffect_interaction_late <- significanteffect_interaction[significanteffect_interaction>50]

pe_coefficients$coef_interaction_pe_valence[significanteffect_interaction_late] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_interaction_pe_valence == min(pe_coefficients$coef_interaction_pe_valence[significanteffect_interaction_late]))
which(pe_coefficients$coef_interaction_pe_valence == min(pe_coefficients$coef_interaction_pe_valence[significanteffect_interaction_late]))*4-704
# minimum at samplepoint 172, i.e. at -16 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_interaction_pe_valence == min(pe_coefficients$coef_interaction_pe_valence[significanteffect_interaction_late]))
                      , c("coef_interaction_pe_valence", "se_interaction_pe_valence", "prob_interaction_pe_valence", "df_interaction_pe_valence","t_interaction_pe_valence", "effectsize_interaction")],3)
round(pe_coefficients[which(pe_coefficients$coef_interaction_pe_valence == min(pe_coefficients$coef_interaction_pe_valence[significanteffect_interaction_late]))
                      , c("coef_interaction_pe_valence", "se_interaction_pe_valence", "prob_interaction_pe_valence", "df_interaction_pe_valence","t_interaction_pe_valence", "effectsize_interaction")],2)

## simple slope of pos PEs around omission onset

significanteffect_pos_pe_late <- significanteffect_pos_pe[significanteffect_pos_pe>50]

pe_coefficients$coef_pos[significanteffect_pos_pe_late] # coefficients are all negatively signed, therefore we look for the minimum
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe_late]))
which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe_late]))*4-704
# minimum at samplepoint 173, i.e. at -12 ms relative to fb omission onset (samplepoint*4-704)

# get relevant info
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe_late]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],3)
round(pe_coefficients[which(pe_coefficients$coef_pos == min(pe_coefficients$coef_pos[significanteffect_pos_pe_late]))
                      , c("coef_pos", "se_pos", "pvalue_pos")],2)







