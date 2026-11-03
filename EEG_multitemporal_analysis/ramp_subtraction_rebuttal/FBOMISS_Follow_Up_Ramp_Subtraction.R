#### FBOmission: Multitemporal Analysis: Follow Up With Ramp Subtraction (Rev3) | Extended Segment (-700 to 1500 ms) ####

### (CW, 11/2026)

##### Routine #####

remove(list = ls()) # clear workspace
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss") # set working directory

#### Read data (Second best BIC (5b/5.1)) ####
data <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_700pre_1500post_anonym_5b.csv") # read concatenated data


#### Exclude participants without variance in the model-derived |PE| ####

# For participants whose fitted learning rates are ~0 (in combination with high
# inverse temperatures), |PE| stays ~0.5 on every trial. Their by-participant 
# PE slopes cannot be estimated, which makes the LMM fits unstable (spikes in 
# the PE coefficient). Criterion: SD of |PE| < 0.05 in omission or display 
# trials (same sample for both conditions).

sd_criterion <- 0.05

pe_sd <- data[!is.na(PE) & omission %in% c(-1, 1),
              .(sd_abs_pe = sd(abs(PE)), n_trials = .N), by = .(id, omission)]

low_var_ids <- unique(pe_sd[sd_abs_pe < sd_criterion, id])

print(pe_sd[id %in% low_var_ids][order(id, omission)]) # SDs of the excluded participants
print(paste0(length(low_var_ids), " participant(s) excluded from the PE analyses: ",
             paste(sort(low_var_ids), collapse = ", ")))

data <- data[!(id %in% low_var_ids)]
print(paste0("remaining participants: ", length(unique(data$id))))

#### Ramp Subtraction ####

## Logic (as suggested by Reviewer 3): for each trial and electrode, fit a
## linear trend to the pre-omission period (-400 to 0 ms), extrapolate it over
## the whole segment and subtract it. If post-onset PE effects merely continue
## the pre-onset ramp, they should disappear in the residuals.
##
## Sample timing: 250 Hz, sample i corresponds to (i * 4 - 704) ms,
## i.e. sample 1 = -700 ms, sample 76 = -400 ms, sample 176 = 0 ms, sample 550 = 1492 ms.

fit_from <- 76   # -400 ms
fit_to   <- 176  #    0 ms
n_samples <- 550

all_electrodes <- c("CP1", "CP2", "P3", "Pz", "P4") # centroparietal cluster

t_all <- (1:n_samples) * 4 - 704          # time in ms for all samples
t_fit <- t_all[fit_from:fit_to]           # time in ms for the fit window
t_fit_c <- t_fit - mean(t_fit)

for (el in all_electrodes) {
  cols <- paste0(el, "_", 1:n_samples)
  Y <- as.matrix(data[, ..cols])                          # trials x samples
  Y_fit <- Y[, fit_from:fit_to, drop = FALSE]
  
  # least-squares line per trial (closed form, vectorised)
  slope <- as.vector(Y_fit %*% t_fit_c) / sum(t_fit_c^2)
  intercept <- rowMeans(Y_fit) - slope * mean(t_fit)
  
  # extrapolated ramp for every sample and residual amplitudes
  ramp <- outer(intercept, rep(1, n_samples)) + outer(slope, t_all)
  data[, (cols) := data.table::as.data.table(Y - ramp)]
}


#### Subset ####
data_presentation <- subset(data, omission==-1)
data_omission <- subset(data, omission==1)

segment_length <- 550 
set.seed(17) # for reproducibility

##### Multitemporal Analysis: Omission Data | Prediction Error | Centroparietal electrodes #####

library(lmerTest) # generates p-values and automatically loads lme4

# subset data
#aggr_data <- data # save full data in other object
data <- data_omission # subset only omission trials
data <- as.data.frame(data) # ensure that data is a data frame (and not a list)

gc()

### Prepare variables

## Prediction Error

# take absolute PE
data$pe_absolute <- abs(data$PE)
hist(data$pe_absolute) #check

# rescale to -1 and 1 for better comparison with categorical variable (valence)
data$pe_absolute <- scales::rescale(data$pe_absolute, to = c(-1, 1))
hist(data$pe_absolute) #check

# center around mean
data$pe_absolute <- data$pe_absolute - mean(data$pe_absolute, na.rm=T)
hist(data$pe_absolute)

# convert to data table (which is more efficient)
data <- data.table::as.data.table(data)

## Valence

hist(data$valence)
# center around mean as suggested by Reviewer 3
data$valence <- data$valence - mean(data$valence, na.rm=T)
hist(data$valence)

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
    
    pe_coefficients[i, "convergence_false"] <- 0
    # enter default value for non-convergence (which is overwritten when a 
    # warning is caught with calling handlers (see below)
    
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
    
    
    
    # compute effectsizes and save (again only of effects with preference level)
    effectsizes <- effectsize::eta_squared(samplepoint_regression)
    
    pe_coefficients[i,"effectsize_valence"] <- effectsizes[which(effectsizes$Parameter=="valence"),"Eta2_partial"]
    pe_coefficients[i,"effectsize_pe"] <- effectsizes[which(effectsizes$Parameter=="pe_absolute"),"Eta2_partial"]
    pe_coefficients[i,"effectsize_interaction"] <-  effectsizes[which(effectsizes$Parameter=="pe_absolute:valence"),"Eta2_partial"]
    
    
    
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
    pe_coefficients[i,"coef_neg"] <- emtrends[which(emtrends$valence<0),c("pe_absolute.trend")]
    pe_coefficients[i,"coef_pos"] <- emtrends[which(emtrends$valence>0),c("pe_absolute.trend")]
    
    pe_coefficients[i,"se_neg"] <- emtrends[which(emtrends$valence<0),c("SE")]
    pe_coefficients[i,"se_pos"] <- emtrends[which(emtrends$valence>0),c("SE")]
    
    pe_coefficients[i,"pvalue_neg"] <- emtrends[which(emtrends$valence<0),c("p.value")]
    pe_coefficients[i,"pvalue_pos"] <- emtrends[which(emtrends$valence>0),c("p.value")]
    
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
  
  pe_coefficients[i, "convergence_false_message"] <<- w$message #  save warning message
  pe_coefficients[i, "convergence_false"] <<- length(w$message) # overwrite 0 for dummy-coded convergence index
  
})


# add number of samplepoints for plotting
pe_coefficients$time <- 1:nrow(pe_coefficients) 

#### Save data in working directory
pe_coefficients_electrode_cluster <- pe_coefficients
data.table::fwrite(pe_coefficients_electrode_cluster, "aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_omission_trials_no_baseline_5b_with_ramp_subtraction.csv", row.names=F)


#### Plot of ramp-corrected results | M5.1 | Omission | Centroparietal ####

pe_coefficients_centroparietal <- pe_coefficients


# General settings (as in the main plot script)
sampling_rate <- 250
sampling_width <- 1000/sampling_rate

segment_start <- 700
segment_end <- 1500
segment_length <- 550

# axis labels
axis_start <- -500
axis_label_distance <- 500

# axis ticks (in samplepoints)
axis_start_sample_points <- (axis_start+segment_start+sampling_width)/sampling_width
axis_distance <- axis_label_distance/sampling_width
zero_line_sample_points <- (segment_start+sampling_width)/sampling_width

ylim <- c(4,-4) # limits of y-axis in all plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/4)
x_lab_text <- 175


#### Open plot device ####
png(filename = paste0("plots/multitemp_results_omission_extended_segment_no_baseline_ramp_corrected_5b_",as.character(Sys.Date()),".png"),
    width=4, height=5, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(2,1), mai=c(.67,.55,.4,0.2), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

#### Plot 1 (top): Effects at Centroparietal Cluster | PE ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction
significanteffect_pe <- which(p.adjust(pe_coefficients$prob_pe, method = "BH") < .05)
significanteffect_valence <- which(p.adjust(pe_coefficients$prob_valence, method = "BH") < .05)
significanteffect_interaction <- which(p.adjust(pe_coefficients$prob_interaction_pe_valence, method = "BH") < .05)

## And only keep effects at 5 consecutive sample points

# custom function
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


## Create index for non-converging and singular-fit models
non_convergence <- which(pe_coefficients$convergence_false > 0)
singular <- which(pe_coefficients$singular > 0)

## Create plot
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "") # empty plot

## Plot significance indicators

# to loop through fixed effects:
effects <- c("significanteffect_pe", "significanteffect_valence", "significanteffect_interaction")

# Define the mixed color
mid_color <- "#4D95A6" # opaque equivalent of #006780 at alpha .7 on white (avoids darkening where significance squares overlap)
col <- c("purple", mid_color, "cornflowerblue")
pch <- c(15, 15, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect)>0){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}

## Add indicators for model convergence and singular fit

# add coral circles for singular fit models
if (length(singular)>0){
  for (i in 1:length(singular)){
    points(singular[i], (ylim[1]-(.0375*abs(diff(ylim)))), pch=19, col="coral", xpd=T)
  }
}

# add golden circles for non-converged model 
if (length(non_convergence)>0){
  for (i in 1:length(non_convergence)){
    points(non_convergence[i], ylim[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
  }
}


## Draw lines for coefficients


lines(pe_coefficients$time, pe_coefficients$coef_interaction_pe_valence, col="cornflowerblue", lwd = 1, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_interaction_pe_valence + pe_coefficients$se_interaction_pe_valence, 
                       rev(pe_coefficients$coef_interaction_pe_valence - pe_coefficients$se_interaction_pe_valence))),
        col=yarrr::transparent("cornflowerblue", trans.val = .7),xpd=T, border=NA)


lines(pe_coefficients$time, pe_coefficients$coef_pe, col="purple", lwd = 1, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                       rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent("purple", trans.val = .7),xpd=T, border=NA)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_valence + pe_coefficients$se_valence, 
                       rev(pe_coefficients$coef_valence - pe_coefficients$se_valence))), col=yarrr::transparent(mid_color, trans.val = .7),xpd=T, border=NA)


lines(pe_coefficients$time, pe_coefficients$coef_valence, col=mid_color, lwd = 1, xpd=T)

## Annotations

# Axes
axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
#text(x_lab_text,(ylim[1] + .3*abs(diff(ylim))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset
# y-axis label
title(ylab = expression(beta~"-coefficient"))

mid_color <- "#4D95A6" # opaque equivalent of #006780 at alpha .7 on white (avoids darkening where significance squares overlap)
legend("topright", c(as.expression(bquote("Absolute PE")),
               as.expression(bquote("Feedback Valence")),
               as.expression(bquote("Absolute PE x Feedback Valence"))), 
       col=c("purple",mid_color, "cornflowerblue"),
       bty="n", xpd=T,
       x.intersp=.7, lty=1, lwd=4, cex=.6)
legend(-20,.5, c("No Convergence", "Singular-fit"),
       col=c("darkgoldenrod1", "coral"), pch=19, bty="n", xpd=T, #title= "Model Diagnostics",
       x.intersp=.7, cex=.7)#, cex=.8) # add legend for line colours


#### Plot 2 (bottom): Simple Slopes of PE at Centroparietal Cluster ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal

## Apply alpha correction
significanteffect_pe <- which(p.adjust(pe_coefficients$prob_pe, method = "BH") < .05)
significanteffect_valence <- which(p.adjust(pe_coefficients$prob_valence, method = "BH") < .05)
significanteffect_interaction <- which(p.adjust(pe_coefficients$prob_interaction_pe_valence, method = "BH") < .05)

## And only keep effects at 5 consecutive sample points

# custom function
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

## Create plot
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "") # empty plot

## Plot significance indicators


# to loop through effect:
effects <- c( "significanteffect_pos_pe",
              "significanteffect_neg_pe", "significanteffect_diff_valence_pe")

# Define the mixed color
mid_color <- "#4D95A6" # opaque equivalent of #006780 at alpha .7 on white (avoids darkening where significance squares overlap)
col <- c("blue", "green3", "cyan1")

pch <- c(15, 15, 15)
#pch <- c(23, 23,23,  19, 19, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect)>0){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.28*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}



## Draw lines for coefficients
lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 1, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 1, xpd=T)
#lines(pe_coefficients$time, pe_coefficients$coeff_diff, col="cyan1", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
                       rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
                       rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)

## Annotations

# Axes
axis(1, at = seq(axis_start_sample_points,segment_length+1,axis_distance), labels = seq(axis_start, segment_end+sampling_width, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim[1] + .4*abs(diff(ylim))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset

# y-axis label
title(ylab = expression(beta~"-coefficient of PE"))

legend(180,-5, c("Positive Outcomes", "Negative Outcomes", "Diff. Test of Simple Slopes"),
       col=c("blue","green3", "cyan1"),bty="n", xpd=T, lty=c(1,1,NA), pch=c(NA,NA,14), lwd=c(4,4,4),#lwd=c(6,6),
       x.intersp=.7, cex=.7) # add legend for line colours

#### Close plot ####
dev.off() # close device
