## FBOmission - Behavioural Data Analysis

## CW, 01/2026

#### General settings & loading of packages ####
remove(list = ls()) # clear workspace
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss")

# for data wrangling
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)

# for statistical tests
library(rstatix)
library(lme4)
library(lmerTest)
library(buildmer)
library(coin)
library(emmeans)

# for plots
library(ggeffects)
library(ggplot2)
library("RColorBrewer")
library(viridis)

#general settings
options(width=500, scipen=6, digits=8)
options("jtools-digits" = 3)
control=glmerControl(optCtrl=list(maxfun=1e6))

set.seed(17)

#### Learning Performance (Choice Accuracy) ####

## Aggregation of behavioural data 

# In the analysis examining accuracy of choices (i.e. evidence of learning),
# only responses for the 70% and 90% stimuli are included as there was no right
# or wrong response for the 50% stimulus. Non-learners are defined as having 
# responded less accurate than chance and will be excluded from all further analyses.

# Read behavioural data
behav <- data.table::fread("aggregated_data/FBOmiss_behaviour_immediate_anonymized.csv")
behav <- as.data.frame(behav)

# Create new columns for feedback valence, type (omission/presentation), learning context (reward/loss), and reward probability

# previous coding of feedback:
# 11: positive presented; 
# -11: negative presented;
# 1: positive omitted;
# -1: negative omitted.

behav$valence <- ifelse(behav$feedback == 11 | behav$feedback == 1, 1, ifelse(behav$feedback == -11 | behav$feedback == -1, -1, NA)) # positive outcome = 1; negative outcome = -1
behav$omission <- ifelse(behav$feedback == 1 | behav$feedback == -1, 1, ifelse(behav$feedback == -11 | behav$feedback == 11, -1, NA)) # omitted feedback = 1; presented feedback = -1
behav$context <- ifelse(behav$feedback == 11 | behav$feedback == -1, 1, ifelse(behav$feedback == -11 | behav$feedback == 1, -1, NA)) # reward context = 1; loss context = -1
behav$rew_prob <- ifelse(behav$stim == 3 | behav$stim == 4, 70, ifelse(behav$stim==1|behav$stim==2, 90, NA))

# Recode "9" to NAs in "correct" and "choice" columns
behav[which(behav$correct == 9), "correct"] <- NA
behav[which(behav$choice == 9), "choice"] <- NA

# Stimulus 5 & 6 had a reward probability of 50%; accordingly, there is no "correct choice".
# Therefore I recode the correct column for these stimuli before evaluating accuracy rates.
behav[which(behav$stim == 5|behav$stim == 6), "correct"] <- NA

# Check whether performance was below chance based on binomial distribution

learnerlist <- data.frame()

# loop through all separate behavioural datafiles:
for (i in 1:length(unique(behav$id))){
  
  # file currently evaluated:
  file <- unique(behav$id)[i]
  
  # total number of trials with responses in current file:
  size <- sum(!is.na(behav$correct[which(behav$id==file)]))
  # number of trials with correct responses in current file:
  n_correct <- sum(behav$correct[which(behav$id==file)], na.rm=T)
  
  # minimum of correct responses (of total number of trials of current data (defined in
  # size)) to not be classified as "false below chance"
  # with a p-value below .05 and a probability to be correct at random of .5:
  n_below_chance <- qbinom(.05, size=size, prob=.5, lower.tail=TRUE)
  
  # fill data frame 
  learnerlist[i,1] <- file
  learnerlist[i,2] <- ifelse(n_correct <= n_below_chance, 1, 0)
  
}


# name columns of new object
names(learnerlist) <- c("id", "below_chance")

# save list of non-learners
non_learner <- learnerlist[which(learnerlist[,2]>0),"id"]

# exclude these participants
if (length(non_learner)>0){
  behav <- behav[-which(behav[,"id"]%in%non_learner),]
}

print(paste0(length(learnerlist[which(learnerlist[,2]>0),"id"]), " participants performed significantly below chance and will be excluded. The resulting samples comprises ", nrow(behav)/480, " participants."))


## Statistical Analysis of accuracy data: Data Preparation

behav$trialno_full <- rep(1:480, times=nrow(behav)/480) #overwrite trialno_full because it starts counting again after half of the experiment

# scale trial number around -1 and 1
behav$trialno_centered <- 2*(behav$trialno_full - min(behav$trialno_full)) / (max(behav$trialno_full) - min(behav$trialno_full))-1

# effect coding for reward probability (only includes learnable stimuli)
behav$rew_prob <- ifelse(behav$rew_prob==70,1,-1)


## GLMM

accuracy_glmer<- glmer(correct ~ trialno_centered*context*rew_prob + (1+trialno_centered*context*rew_prob|id), family=binomial, data=behav, control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=1e5)))

summary(accuracy_glmer)
performance::r2(accuracy_glmer,  tolerance = 1e-10) # model fit/explanatory power


#### Resolve interaction

# Plot interaction (with odds logit)
emmeans::emmip(accuracy_glmer,rew_prob ~ trialno_centered|  context, CIs = TRUE, cov.reduce = range, lmer.df="asymp")

# Resolve interaction
emmeans::emtrends(accuracy_glmer, 
                  specs=pairwise~rew_prob+context| trialno_centered, 
                  var="trialno_centered",
                  infer=TRUE,
                  adj="bonf", lmer.df="asymp")

# All simple slopes are significant. However, contrasts between slopes show that
# the slope for trial number is steeper for the 90% stimulus in the Get Reward 
# Context vs the 90% stimulus in Avoid Loss Context (p=.0083), vs 70% stimulus 
# in Avoid Loss Context (p=.0126), but not vs the 70% stimulus in Get Reward 
# Context (p = .1521).

#### Visualization

# create vector for range of predicted values
trial_seq <- seq(-1, 1, length.out = 480)

# create predictions based on the glmm on reponse scale
emm <- emmeans(accuracy_glmer,
               ~ trialno_centered | context * rew_prob,
               at = list(trialno_centered = trial_seq),
               type = "response")

# convert to dataframe
plotdata <- as.data.frame(emm)

# transform to percent
plotdata$accuracy_percent <- plotdata$prob * 100
plotdata$lower_ci <- plotdata$asymp.LCL * 100
plotdata$upper_ci <- plotdata$asymp.UCL * 100

# create trial bins for aggregation of raw data
n_bins <- 8
behav$trial_bin <- cut(behav$trialno_centered,
                       breaks = seq(-1, 1, length.out = n_bins + 1),
                       include.lowest = TRUE)

# create empty dataframe
agg_data <- data.frame()


for (ctx in c(-1, 1)) { # loop through contexts 
  for (rw in c(-1, 1)) { # and outcome probabilities
    
    # subset
    subset_data <- behav[behav$context == ctx & behav$rew_prob == rw, ]
    # create vector of bins for loop
    bins <- unique(subset_data[!is.na(subset_data$trial_bin),"trial_bin"])
    
    for (bin in bins) { # loop through bins
      
      # subset
      bin_data <- subset_data[subset_data$trial_bin == bin, c("id", "correct")]
      
      # aggregate at first per participant
      means_subject <- aggregate(correct ~ id, bin_data, mean)
      # save average column as vector
      bin_data <- means_subject$correct
      n <- length(bin_data)
      
      if (n > 1) {
        # average over participants
        m <- mean(bin_data, na.rm=T)
        # calculate standard error
        se <- sd(bin_data, na.rm=T) / sqrt(n)
        # and ci boundaries
        lower <- (m - 1.96 * se) * 100
        upper <- (m + 1.96 * se) * 100
        m <- m * 100
        
        # calculate bin center for x position (jittered it now...)
        bin_center <- as.numeric(sub(".*,", "", gsub("[^,0-9.-]", "", bin))) - 1 / n_bins
        
        # bind to agg_data
        agg_data <- rbind(agg_data, data.frame(trial_bin = bin, rew_prob = rw, 
                                               context = ctx, correct_mean = m, lower = lower, upper =
                                                 upper, trialno_centered = bin_center))
      }
    }
  }
}


# create plot
png("plots/predicted_accuracy.png", width = 3.8, height = 1.9, units = "in", res = 2400)

par(mfrow = c(1, 2), mai = c(0.38, 0.4, 0.2, 0.05), mgp = c(1.4, 0.35, 0), tcl = -0.25, cex = 0.8)

for (ctx in c(1, -1)) { # loop through contexts
  
  # subset predicted and observed data
  plot_data <- plotdata[plotdata$context == ctx, ]
  obs_data <- agg_data[agg_data$context == ctx, ]
  
  ctx_label <- ifelse(ctx == 1, "Get Reward", "Avoid Loss")
  
  # create empty plot 
  if (ctx==1){ #(for first plot with y-axis label)
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "Choice Accuracy (%)",
         xaxt = "n", yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
  } else{
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "", xaxt = "n",
         yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
  }
  
  
  # plot axes
  x_labels_original <- c(1, 120, 240, 360, 480) # prepare x-axis labels 
  x_vals_centered <- (x_labels_original - 240.5) / 239.5 # mapped to centered trialvariable
  axis(1, at = x_vals_centered, labels = x_labels_original, cex.axis = 0.8, mgp=c(1, 0.1, 0))
  text(0,-22, "Trial", xpd=T, cex=.8)
  axis(2, las=2, cex.axis = 0.8)
  
  # define colors for outcome probabilities
  col_90 <- "#1B9E77"   # Teal
  col_70 <- "#6C7A89"   # Slate gray
  
  for (rw in c(-1, 1)) { # loop through outcome probabilities
    
    # assign color and linetype (green dashed line for 90% stim, grey solid line for 70% stim)
    col <- ifelse(rw == -1, col_90, col_70)
    lty <- ifelse(rw == -1, 2, 1)
    
    # subset glmm-predicted and observed data
    pd <- plot_data[plot_data$rew_prob == rw, ]
    od <- obs_data[obs_data$rew_prob == rw, ]
    
    # plot predicted data
    
    # shaded areas for confidence intervals
    polygon(c(pd$trialno_centered, rev(pd$trialno_centered)),
            c(pd$lower_ci, rev(pd$upper_ci)), col = adjustcolor(col, alpha.f = 0.2), border = NA)
    # lines for predicted accuracy
    lines(pd$trialno_centered, pd$accuracy_percent, col = col, lty = lty, lwd = 1.8)
    
    # plot observed data
    
    # prepare jitter
    jitter <- ifelse(rw==-1, .035, -.035)
    od$trialno_centered <- od$trialno_centered+jitter
    # error bars for confidence intervals
    arrows(od$trialno_centered, od$lower, od$trialno_centered, od$upper,
           angle = 90, code = 3, length = 0.02, col = col, lwd = 0.8)
    # averaged accuracy as points
    points(od$trialno_centered, od$correct_mean, col = col, pch = 16, cex = 0.5)
    
    
  }
  
  # add legend in second plot only
  if (ctx==-1) {legend("bottomright", legend = c("90%", "70%"),
                       col = c(col_90, col_70), lty = c(2, 1), lwd = 1.5, #pch = 16,
                       bty = "n", cex = 0.60)
    text(-0.31,8, expression(P("Pos. Outcome" * "|" * "Correct")), xpd=T, cex=.65)}  
}

dev.off()


#### Analysis of accuracy data from simulated_data ####

behav_sim <- read.csv("behavioural_data/immediate/comp_model_fit_export/FBOmiss_behav_recovered.csv", sep=",", header=T)
# no filenames etc. in the file (no need to anonymize)

behav_obs <- behav # rename observed data
behav <- behav_sim

# recode variables such that I can use the same plot script as above for observed data

# recode and create variables for valence, omission, context, and reward probability
behav$valence <- ifelse(behav$feedback == 1 , 1, ifelse(behav$feedback == 0, -1, NA)) # positive outcome = 1; negative outcome = -1
behav$omission <- ifelse(behav$feedback_type == 1, -1, ifelse(behav$feedback_type == 0, 1, NA)) # omitted feedback = 1; presented feedback = -1
behav$context <- ifelse(behav$feedback == 1 & behav$feedback_type == 1 | behav$feedback == 0 & behav$feedback_type == 0 , 1, ifelse(behav$feedback == 0 & behav$feedback_type == 1 | behav$feedback == 1 & behav$feedback_type == 0, -1, NA)) # reward context = 1; loss context = -1
behav$rew_prob <- ifelse(behav$stim == 3 | behav$stim == 4, 70, ifelse(behav$stim==1|behav$stim==2, 90, NA))

# rename accuracy variable like in observed data
names(behav)[which(names(behav)=="accuracy")] <- "correct"

# Stimulus 5 & 6 had a reward probability of 50%; accordingly, there is no "correct choice".
# Therefore I recode the correct column for these stimuli before evaluating accuracy rates.
behav[which(behav$stim == 5|behav$stim == 6), "correct"] <- NA

# scale trial number around -1 and 1
behav$trialno_centered <- 2*(behav$trial - min(behav$trial)) / (max(behav$trial) - min(behav$trial))-1

# effect coding for reward probability (only includes learnable stimuli)
behav$rew_prob <- ifelse(behav$rew_prob==70,1,-1)


## GLMM

accuracy_glmer_recovered <- lme4::glmer(correct ~ trialno_centered*context*rew_prob + (1+trialno_centered*context*rew_prob|subj:sim_iteration), family=binomial, data=behav, control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=1e5)))
summary(accuracy_glmer_recovered)
performance::r2(accuracy_glmer_recovered, tolerance = 1e-10) # model fit/explanatory power

#### Resolve interaction

# Plot interaction (with odds logit)
emmeans::emmip(accuracy_glmer_recovered,rew_prob ~ trialno_centered|  context, CIs = TRUE, cov.reduce = range, lmer.df="asymp")

# Resolve interaction
emmeans::emtrends(accuracy_glmer_recovered, 
                  specs=pairwise~rew_prob+context| trialno_centered, 
                  var="trialno_centered",
                  infer=TRUE,
                  adj="bonf", lmer.df="asymp")

#### Visualization

# create vector for range of predicted values
trial_seq <- seq(-1, 1, length.out = 480)

# create predictions based on the glmm on reponse scale
emm <- emmeans(accuracy_glmer_recovered,
               ~ trialno_centered | context * rew_prob,
               at = list(trialno_centered = trial_seq),
               type = "response")

# convert to dataframe
plotdata <- as.data.frame(emm)

# transform to percent
plotdata$accuracy_percent <- plotdata$prob * 100
plotdata$lower_ci <- plotdata$asymp.LCL * 100
plotdata$upper_ci <- plotdata$asymp.UCL * 100

# create trial bins for aggregation of raw data
n_bins <- 8
behav$trial_bin <- cut(behav$trialno_centered,
                       breaks = seq(-1, 1, length.out = n_bins + 1),
                       include.lowest = TRUE)

# create empty dataframe
agg_data <- data.frame()


for (ctx in c(-1, 1)) { # loop through contexts 
  for (rw in c(-1, 1)) { # and outcome probabilities
    
    # subset
    subset_data <- behav[behav$context == ctx & behav$rew_prob == rw, ]
    # create vector of bins for loop
    bins <- unique(subset_data[!is.na(subset_data$trial_bin),"trial_bin"])
    
    for (bin in bins) { # loop through bins
      
      # subset
      bin_data <- subset_data[subset_data$trial_bin == bin, c("subj","sim_iteration", "correct")]
      
      # aggregate at first across simulations per participant
      means_sims <- aggregate(correct ~ subj + sim_iteration, bin_data, mean)
      
      # then per participant
      means_subject <- aggregate(correct ~ subj, means_sims, mean)
      
      # save average column as vector
      bin_data <- means_subject$correct
      n <- length(bin_data)
      
      if (n > 1) {
        # average over participants
        m <- mean(bin_data, na.rm=T)
        # calculate standard error
        se <- sd(bin_data, na.rm=T) / sqrt(n)
        # and ci boundaries
        lower <- (m - 1.96 * se) * 100
        upper <- (m + 1.96 * se) * 100
        m <- m * 100
        
        # calculate bin center for x position (jittered it now...)
        bin_center <- as.numeric(sub(".*,", "", gsub("[^,0-9.-]", "", bin))) - 1 / n_bins
        
        # bind to agg_data
        agg_data <- rbind(agg_data, data.frame(trial_bin = bin, rew_prob = rw, 
                                               context = ctx, correct_mean = m, lower = lower, upper =
                                                 upper, trialno_centered = bin_center))
      }
    }
  }
}


# create plot
png("plots/accuracy_posterior_prediction.png", width = 3.8, height = 1.9, units = "in", res = 2400)

par(mfrow = c(1, 2), mai = c(0.38, 0.4, 0.2, 0.05), mgp = c(1.4, 0.35, 0), tcl = -0.25, cex = 0.8)

for (ctx in c(1, -1)) { # loop through contexts
  
  # subset predicted and observed data
  plot_data <- plotdata[plotdata$context == ctx, ]
  obs_data <- agg_data[agg_data$context == ctx, ]
  
  ctx_label <- ifelse(ctx == 1, "Get Reward", "Avoid Loss")
  
  # create empty plot (for first plot with y-axis label)
  if (ctx==1){
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "Choice Accuracy (%)",
         xaxt = "n", yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8)
  } else{
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "", xaxt = "n",
         yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8)
    # add legend
    legend("bottomright", legend = c("90%", "70%"),
           col = c(col_90, col_70), lty = c(2, 1), lwd = 1.5, #pch = 16,
           bty = "n", cex = 0.60)
    text(-0.31,8, expression(P("Pos. Outcome" * "|" * "Correct")), xpd=T, cex=.65)
  }
  
  # plot axes
  x_labels_original <- c(1, 120, 240, 360, 480) # prepare x-axis labels 
  x_vals_centered <- (x_labels_original - 240.5) / 239.5 # mapped to centered trialvariable
  axis(1, at = x_vals_centered, labels = x_labels_original, cex.axis = 0.8, mgp=c(1, 0.1, 0))
  text(0,-22, "Trial", xpd=T, cex=.8)
  axis(2, las=2, cex.axis = 0.8)
  
  # define colors for outcome probabilities
  col_90 <- "#1B9E77"   # Teal
  col_70 <- "#6C7A89"   # Slate gray
  
  for (rw in c(-1, 1)) { # loop through outcome probabilities
    col <- ifelse(rw == -1, col_90, col_70)
    lty <- ifelse(rw == -1, 2, 1)
    
    pd <- plot_data[plot_data$rew_prob == rw, ]
    od <- obs_data[obs_data$rew_prob == rw, ]
    
    # plot redicted data
    
    # shaded areas for confidence intervals
    polygon(c(pd$trialno_centered, rev(pd$trialno_centered)),
            c(pd$lower_ci, rev(pd$upper_ci)), col = adjustcolor(col, alpha.f = 0.2), border = NA)
    # lines for predicted accuracy
    lines(pd$trialno_centered, pd$accuracy_percent, col = col, lty = lty, lwd = 1.8)
    
    # plot observed data
    
    # prepare jitter
    jitter <- ifelse(rw==-1, .035, -.035)
    od$trialno_centered <- od$trialno_centered+jitter
    # error bars for confidence intervals
    arrows(od$trialno_centered, od$lower, od$trialno_centered, od$upper,
           angle = 90, code = 3, length = 0.02, col = col, lwd = 0.8)
    # averaged accuracy as points
    points(od$trialno_centered, od$correct_mean, col = col, pch = 16, cex = 0.5)
    
  }
}

dev.off()

#### Analysis of (Computationally-modelled) Learning Rates ####

## Data preparation

# read data
parameter <- read.csv("behavioural_data/immediate/comp_model_fit_export/FBOmiss_learning_parameter.csv")

# show descriptive statistics for parameters of rl model
params <- names(parameter)[3:15]
psych::describeBy(parameter[,params])

## reduce data (to learning rates, and participant ids):
parameter <- parameter[,c(3:14,16)]

# and reshape
data_learning_rates <- reshape2::melt(parameter, id = c("id"))

# recode variables
data_learning_rates$learning_rate <- as.numeric(data_learning_rates$value)

data_learning_rates$rew_prob <- substr(data_learning_rates$variable, 7,8)
data_learning_rates[grep("pos", data_learning_rates$variable), "valence"] <- "pos"
data_learning_rates[grep("neg", data_learning_rates$variable), "valence"] <- "neg"

data_learning_rates[grep("omitted", data_learning_rates$variable), "type"] <- "omission"
data_learning_rates[grep("presented", data_learning_rates$variable), "type"] <- "presentation"

## Transform data for approximating normal distribution

### Original distribution

# Plot distribution
qqnorm(data_learning_rates$learning_rate, pch = 1, frame = FALSE)
qqline(data_learning_rates$learning_rate, col = "steelblue", lwd = 2)

hist(data_learning_rates$learning_rate)
plot(density(data_learning_rates$learning_rate))

### Log transformed distribution

data_learning_rates$learning_log <- log(data_learning_rates$learning_rate)

# Plot distribution
qqnorm(data_learning_rates$learning_log, pch = 1, frame = FALSE)
qqline(data_learning_rates$learning_log, col = "steelblue", lwd = 2)

hist(data_learning_rates$learning_log)
plot(density(data_learning_rates$learning_log))

### Square root transformed distribution

data_learning_rates$learning_square_root <- data_learning_rates$learning_rate^(1/2)

# Plot distribution
qqnorm(data_learning_rates$learning_square_root, pch = 1, frame = FALSE)
qqline(data_learning_rates$learning_square_root, col = "steelblue", lwd = 2)

hist(data_learning_rates$learning_square_root)
plot(density(data_learning_rates$learning_square_root))

### Cube-root transformed distribution

data_learning_rates$learning_rate_cube_root <- data_learning_rates$learning_rate^(1/3)

# Plot distribution
qqnorm(data_learning_rates$learning_rate_cube_root, pch = 1, frame = FALSE)
qqline(data_learning_rates$learning_rate_cube_root, col = "steelblue", lwd = 2)

hist(data_learning_rates$learning_rate_cube_root)
plot(density(data_learning_rates$learning_rate_cube_root))


# Cube-root transformation is chosen based on the normal qqplots, therefore I replace the variable in the dataset:
data_learning_rates$learning_original <- data_learning_rates$learning_rate
data_learning_rates$learning_rate <- data_learning_rates$learning_rate^(1/3)

## Statistical analysis with linear mixed effects models: Data Preparation

data <- data_learning_rates

# recode variables
data$valence <- ifelse(data$valence == "pos", 1, ifelse(data$valence == "neg", -1, NA))
data$omission <- ifelse(data$type == "omission", 1, ifelse(data$type == "presentation", -1, NA))
data$rew_prob <- ifelse(data$rew_prob == "90", 1, ifelse(data$rew_prob == "70", -1, NA)) # only 70% and 90% included as done for accuracy

data$context <- ifelse((data$valence==1 & data$omission==-1) | (data$valence==-1 & data$omission==1), 1,
                       ifelse((data$valence==1 & data$omission==1) | (data$valence==-1 & data$omission==-1), -1,NA))
# reward context = 1; loss context = -1

## LME

model_learning_rates <- lmer(learning_rate ~ rew_prob*valence*omission + 
                               (1|id),
                             data=data, REML=T, 
                             control=lmerControl(optimizer='bobyqa', optCtr = list(maxfun = 1e9)))
summary(model_learning_rates)
performance::r2(model_learning_rates,  tolerance = 1e-10) 

# inspect distribution of residuals
qqnorm(residuals(model_learning_rates))
qqline(residuals(model_learning_rates))
# looks ok

### Visualize

library(ggeffects)
library(ggplot2)

plotdata <- ggeffect(model_learning_rates, terms = c("omission", "valence", "rew_prob"))

# Rename columns
names(plotdata)[1] <- "omission"
names(plotdata)[6] <- "valence"
names(plotdata)[7] <- "rew_prob"

# Define facetvariable + factorize levels
plotdata$facet_group <- ifelse(plotdata$omission==-1,"Display","Omission")
plotdata$rew_prob <- factor(plotdata$rew_prob, levels = c(-1, 1), labels = c("70%", "90%"))
plotdata$valence <- factor(plotdata$valence, levels = c(-1, 1), labels = c("Negative", "Positive"))

# Split for facetting
facet_split <- split(plotdata, plotdata$facet_group)

# define colors
valence_colors <- c("blue", "green3")
names(valence_colors) <- c("Positive", "Negative")

# Plotting order
facet_order <- c("Display", "Omission")

png("plots/learning_rates_predicted.png", width = 3.8, height = 1.9, units = "in", res = 2400)

par(mfrow = c(1, 2), mai = c(0.38, 0.4, 0.3, 0.05), mgp = c(1.6, 0.35, 0), tcl = -0.25, cex = 0.8)

for (facet_name in facet_order) {
  df <- facet_split[[facet_name]]
  
  # Subset raw data for this facet
  raw_subset <- subset(data, 
                       (valence == 1 & omission == -1 & facet_name == "Display") |
                         (valence == -1 & omission == 1 & facet_name == "Omission") |
                         (valence == -1 & omission == -1 & facet_name == "Display") |
                         (valence == 1 & omission == 1 & facet_name == "Omission")
  )
  
  # Map categorical vars for consistency
  raw_subset$valence <- factor(raw_subset$valence, levels = c(-1, 1), labels = c("Negative", "Positive"))
  raw_subset$omission <- ifelse(raw_subset$omission == -1, "Displayed", "Omitted")
  raw_subset$rew_prob <- factor(raw_subset$rew_prob, levels = c(-1, 1), labels = c("70%", "90%"))
  raw_subset$x <- as.numeric(raw_subset$rew_prob) +
    ifelse(raw_subset$valence == "Positive", -0.15, 0.15) +
    runif(nrow(raw_subset), -0.05, 0.05)  # jitter
  
  # Base plot
  if (facet_name == "Display"){
    # Create empty plot
    plot(1, type = "n", ylim = c(0, 1), xlim = c(0.5, 2.5),
         xaxt = "n", xlab = "", ylab = "Learning Rate",
         main = facet_name, las = 1, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
  } else{ # omit y-lab annotation for second plot
    plot(1, type = "n", ylim = c(0, 1), xlim = c(0.5, 2.5),
         xaxt = "n", xlab = "", ylab = "",
         main = facet_name, las = 1, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
  }
  # Raw points
  points(raw_subset$x, raw_subset$learning_original,
         pch = 16,
         col = adjustcolor(valence_colors[as.character(raw_subset$valence)], alpha.f = 0.4),#previously .05
         cex = 0.5)
  
  # Error bars
  x_vals <- as.numeric(df$rew_prob)
  offset <- ifelse(df$valence == "Positive", -0.15, 0.15)
  x_dodge <- x_vals + offset
  
  arrows(x0 = x_dodge, y0 = df$conf.low, x1 = x_dodge, y1 = df$conf.high,
         angle = 90, code = 3, length = 0.05,
         col = valence_colors[as.character(df$valence)])
  
  # Predicted points
  points(x_dodge, df$predicted,
         pch = 16,
         col = valence_colors[as.character(df$valence)],
         cex = 0.8)
  
  # Axis & label
  axis(1, at = 1:2, labels = levels(df$rew_prob))
  text(1.5, -0.31, "P(Pos. Outcome|Correct)", xpd = TRUE, cex = .9)
  
  
}

dev.off()

# plot legend separately (will be put together outside of R)
png("plots/learning_rates_predicted_legend.png", width = 1.9, height = 1.9, units = "in", res = 2400)

par(mfrow = c(1, 1), mai = c(0.38, 0.4, 0.3, 0.05), mgp = c(1.6, 0.35, 0), tcl = -0.25, cex = 0.8)

plot(1, type = "n", ylim = c(0, 1), xlim = c(0.5, 2.5),
     axes=F,xaxt = "n", xlab = "", ylab = "",
     las = 1, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")

legend("topright", legend = c("Positive", "Negative"),
       pch = 16, col = valence_colors, pt.cex = 0.9,
       bty = "n", cex = .8)


### Resolve interactions

# Simple effects
valence_prob_omission <- emmeans::emmeans(model_learning_rates, 
                                          specs=  ~ valence+rew_prob|omission, 
                                          var="valence",
                                          infer=TRUE,
                                          adj="none", lmer.df="asymp")

pairs(valence_prob_omission)
pairs(pairs(valence_prob_omission), by=NULL)


#### Comparison of fitted with recovered parameters ####

# load fitted (anonymized) parameter
parameter_fitted <- read.csv("behavioural_data/immediate/comp_model_fit_export/FBOmiss_learning_parameter.csv")

# load recovered parameter
parameter_recovered <- read.csv("behavioural_data/immediate/comp_model_fit_export/FBOmiss_learning_parameter_recovered.csv")
names(parameter_recovered)[names(parameter_recovered)=="filename"] <- "subj"

# create list with free parameter
list_parameter <- names(parameter_fitted)
list_parameter <- list_parameter[-c(which(list_parameter=="filename"), which(list_parameter=="BIC"), which(list_parameter=="subj"))]

correlations <- data.frame()
correlations[1:25,1] <- 1:25
names(correlations)[1] <- "sim"

for (j in 1:length(list_parameter)) {
  
  # container to store all 25 recovered vectors
  recovered_mat <- matrix(NA, nrow = nrow(parameter_recovered), ncol = 25)
  
  # loop through simulations
  for (i in 1:25) {
    
    colname <- ifelse(i < 10,
                      paste0(list_parameter[j], "_0", i),
                      paste0(list_parameter[j], "_", i))
    
    if (list_parameter[j] != "beta") { # transform learning rates
      recovered_mat[, i] <- parameter_recovered[, colname]^(1/3)
    } else {
      recovered_mat[, i] <- parameter_recovered[, colname]
    }
  }
  
  # average values across simulations
  recovered_avg <- rowMeans(recovered_mat)
  
  # get corresponding fitted values
  
  if (list_parameter[j] != "beta") {
    fitted <- parameter_fitted[, list_parameter[j]]^(1/3)
    xlim_use <- c(0, 1)
    ylim_use <- c(0, 1)
    
    valence <- ifelse(grepl("neg",list_parameter[j]), "Negative", "Positive")
    rew_prob <- ifelse(grepl("50",list_parameter[j]), "50 %", ifelse(grepl("70",list_parameter[j]),"70 %", "90%"))
    appearance <- ifelse(grepl("omitted",list_parameter[j]), "Omission", "Display")
    if (!exists("plot_title")){
      plot_title <- paste0("LR (", valence, ", ", appearance, ", ", rew_prob, ")")
    }
    
  } else {
    fitted <- parameter_fitted[, list_parameter[j]]
    xlim_use <- c(0, 100)
    ylim_use <- c(0, 100)
    if (!exists("plot_title")){
      plot_title <- "Beta"
    }
  }
  
  # store correlation to write in plot title
  corr_coefficient <- cor(recovered_avg, fitted)
  correlations[j, 2] <- corr_coefficient
  names(correlations)[2] <- list_parameter[j]
  
  # create plot
  png_filename <- file.path(paste0("plots/scatter_", list_parameter[j], ".png"))  
  png(png_filename, width = 1.9, height = 1.9, units = "in", res = 2400)
  
  par(mfrow = c(1, 1), mai = c(0.3, 0.35, 0.3, 0.05), mgp = c(1.6, 0.3, 0), tcl = -0.15, cex = 0.8)
  
  plot(
    recovered_avg, fitted,
    xlab = "",  # averaged across 25 simulations
    ylab = "",
    main = paste0(plot_title, "\n r = ", round(corr_coefficient, 2)),
    #main = "", #paste0(plot_title, "\n r = ", round(corr_coefficient, 2)),
    #main = bquote(alpha ~ .(plot_title) ~ "; r = " ~ .(round(corr_coefficient, 3))),
    cex.main=.8,
    xlim = xlim_use,
    ylim = ylim_use,
    pch = 19,
    col = rgb(0, 0, 0, 0.5),
    cex = 1.1, lwd=0,
    cex.axis = .7, las=1 , xaxt="n" #, ann=F
  )
  
  
  title(ylab="Fitted", line=1.3)
  title(xlab="Simulated", line=1.0)
  
  if (list_parameter[j] != "beta") {
    
    axis(1, at=seq(0,1, 0.2), las=1, cex.axis=.7, mgp=c(1.6,.001,0))
    
  }else{ # adjust axis for beta
    
    axis(1, at=seq(0,100, 20), las=1, cex.axis=.7, mgp=c(1.6,.001,0))
    
  }
  
  
  # add regression line
  abline(lm(fitted ~ recovered_avg), col = "darkgreen", lwd = 2)
  # add line for perfect correlation
  abline(a = 0, b = 1, col = "darkgray", lwd = 2, lty = 2)
  
  dev.off()
  
  print(plot_title)
  rm(plot_title)
  
}


#### Visualize Action Value Estimates ####

# Read action values derived from reinforcement learning models based on the empirical data
pe_data <- read.csv("behavioural_data/immediate/comp_model_fit_export/FBOmiss_immediate_Q_values_and_PEs.csv")
# Rename filename column
pe_data$id <- pe_data$filename

# get observed behavioral data
behav <- behav_obs

# Add trial number variable in both data sets (behavioural data and simulated data)
behav$trial <- rep(1:480, times=(nrow(behav)/480))
pe_data$trial <- rep(1:480, times=(nrow(pe_data)/480))

# Merge data by id and trial number
data <- plyr::join(pe_data, behav, by = c("id", "trial"))

# Add variable for reward probability
data$rew_prob <- ifelse(data$stim == 3 | data$stim == 4, 70, ifelse(data$stim==1|data$stim==2, 90, ifelse(data$stim==5 | data$stim==6, 50, NA)))

# Save data with another for later reuse
data_all <- data

## Add columns in which the action value estimates of the correct responses are transferred (separately for reward probabilities and learning context)

data_sorted <- data.frame()

for (i in 1:length(unique(data$id))){
  
  # subset current subject in loop
  currvp <- unique(data$id)[i]
  currdata <- subset(data, id == currvp)
  
  # create new columns
  newcolumns <- data.frame(matrix(0, ncol = 12, nrow = nrow(currdata)))
  names(newcolumns) <- c("Q_90_reward_corr", "Q_90_loss_corr","Q_70_reward_corr","Q_70_loss_corr", "Q_50_reward_corr", "Q_50_loss_corr",
                         "Q_90_reward_incorr", "Q_90_loss_incorr","Q_70_reward_incorr","Q_70_loss_incorr", "Q_50_reward_incorr", "Q_50_loss_incorr")
  
  # append to current temp data
  currdata <- cbind(currdata, newcolumns)
  
  # Define number of columns of stimuli to sort
  num <- 6
  
  # Loop through each stimulus
  for (j in 1:num) {
    
    # Find the first occurrence where stim is current stim in loop
    first_index <- which(currdata$stim == j & !is.na(currdata$choice))[1]
    
    # Get the values for choice and correct at that index
    choice <- currdata[first_index, "choice"]
    correct <- currdata[first_index, "correct"]
    context <- currdata[first_index, "context"]
    
    
    # Construct the column names dynamically
    
    if (j==1|j==2){
      
      if (context == 1){
        Q_corr <- paste0("Q_90_reward_corr")
        Q_incorr <- paste0("Q_90_reward_incorr")
      } else if (context == -1){
        Q_corr <- paste0("Q_90_loss_corr")
        Q_incorr <- paste0("Q_90_loss_incorr")
      }
    } else if (j==3|j==4){
      
      if (context == 1){
        Q_corr <- paste0("Q_70_reward_corr")
        Q_incorr <- paste0("Q_70_reward_incorr")
      } else if (context == -1){
        Q_corr <- paste0("Q_70_loss_corr")
        Q_incorr <- paste0("Q_70_loss_incorr")
      }
    } else if (j==5|j==6){
      
      if (context == 1){
        Q_corr <- paste0("Q_50_reward_corr")
        Q_incorr <- paste0("Q_50_reward_incorr")
      } else if (context == -1){
        Q_corr <- paste0("Q_50_loss_corr")
        Q_incorr <- paste0("Q_50_loss_incorr")
      }
    }
    
    Q_l <- paste0("Q", j, "l")
    Q_r <- paste0("Q", j, "r")
    
    # Determine the value to copy based on the conditions
    if (is.na(choice == 1 & correct == 1)| is.na(choice == 2 & correct == 1) |
        is.na(choice == 1 & correct == 0)| is.na(choice == 2 & correct == 0)){
      currdata[[Q_corr]] <- NA
    } else if (choice == 1 & correct == 1) {
      currdata[[Q_corr]] <- currdata[[Q_l]]
      currdata[[Q_incorr]] <- currdata[[Q_r]]
    } else if (choice == 1 & correct == 0) {
      currdata[[Q_corr]] <- currdata[[Q_r]]
      currdata[[Q_incorr]] <- currdata[[Q_l]]
    } else if (choice == 2 & correct == 1) {
      currdata[[Q_corr]] <- currdata[[Q_r]]
      currdata[[Q_incorr]] <- currdata[[Q_l]]
    } else if (choice == 2 & correct == 0) {
      currdata[[Q_corr]] <- currdata[[Q_l]]
      currdata[[Q_incorr]] <- currdata[[Q_r]]
    }
    if (j==5|j==6){
      # for the 50% stimuli, I choose right button choice as correct (as there is no "real" right or wrong)
      currdata[[Q_corr]] <- currdata[[Q_r]]
      currdata[[Q_incorr]] <- currdata[[Q_l]]
    }
    
  }
  
  data_sorted <- rbind(data_sorted, currdata)
  
}

data <- data_sorted
data.table::fwrite(data,"aggregated_data/behav_and_pe_data_concatenated.csv", row.names=F) # save data in folder
pe_data <- data # rename in the environment for later merging of scripts with eeg_data

#### Set up axes & define function for expectationplots 

# Set up axes

# save number of possible total trials of one data set for x-axis:
ntrials <- nrow(subset(data, filename==unique(data$filename)[1]))
xaxis <- seq(1, ntrials) # x-axis
yaxis <- c(0,1) # y-axis (here: for expectation/action values ranging between 0 and zero)

# Define customized plot function
expectationplot <- function(x, ylim=yaxis, cols=viridis_pal()(ncol(x)),
                            xlab="Trial", ylab=NA, ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(xaxis), ylim=yaxis, axes=FALSE,
       xlab=xlab, ylab=ylab)
  
  #axis(1) # x-axis
  
  axis(1, at=seq(0, nrow(x), 80)) # x-axis
  #axis(1, at=0:nrow(x)) # x-axis
  
  for (c in 1:ncol(x)){
    lines(xaxis, x[,c], col=cols[c], lwd=2)
  }
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(yaxis[1],yaxis[2], 0.2), pos=0, las=2)
  ## small intermediate axis ticks
  axis(2, at=0:1, pos=0, labels=NA, tcl=par("tcl")/2)
  
  
}

##### Subsample data of conditions and average Qs

columns <- c("Q_90_reward_corr", "Q_90_loss_corr","Q_70_reward_corr","Q_70_loss_corr", "Q_50_reward_corr", "Q_50_loss_corr",
             "Q_90_reward_incorr", "Q_90_loss_incorr","Q_70_reward_incorr","Q_70_loss_incorr", "Q_50_reward_incorr", "Q_50_loss_incorr")

# Create empty data frame with 5 columns and 200 rows
data_q <- data.frame(matrix(0, ncol = length(columns), nrow = ntrials)) # ncol = number of different Q-values, nrow = number of trials

# Loop through the five stimuli and create averages of reward expectations
for (i in 1:length(columns)){
  
  stim <- columns[i]
  
  data_q_temp <- data[, stim] # subset column for stimulus
  data_q_temp <- as.data.frame(split(data_q_temp, 1:ntrials)) # split in trials
  # columns: trials, rows: single subjects
  data_q_temp <- colMeans(data_q_temp, na.rm=T) #average
  data_q[,i] <- data_q_temp # insert in empty data frame column i
  colnames(data_q)[i] <- stim # name column
  
  
}


#### Create plots

# open a plot device
png("plots/action_values_correct.png", width=7, height=3, unit="in", res=2400) #large enough

par(mfcol=c(1,2), mai=c(.7,.7,.3,.05), mgp=c(1.3,.4,0), tcl=-.25)

# Get Reward Context

plot_subset <- data_q[,c(which(names(data_q)=="Q_90_reward_corr"),which(names(data_q)=="Q_70_reward_corr"),which(names(data_q)=="Q_50_reward_corr"))]

expectationplot(plot_subset, ylab="action value estimates")

#add also the values of the incorrect choices
plot_subset <- data_q[,c(which(names(data_q)=="Q_90_reward_incorr"),which(names(data_q)=="Q_70_reward_incorr"),which(names(data_q)=="Q_50_reward_incorr"))]

for (c in 1:ncol(plot_subset)){
  lines(xaxis, plot_subset[,c], col=viridis_pal()(ncol(plot_subset))[c], lwd=2, lty=2) 
}

# insert dashed lines for programmed reward-probability
lines(x=c(0,480), y=c(0.5,0.5), col=viridis_pal()(3)[3], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.7,0.7), col=viridis_pal()(3)[2], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.9,0.9), col=viridis_pal()(3)[1], lty=2, lwd=1)

text(240, 1.1, "Get Reward", xpd=T, cex=1.2)

# part of the legend
text(345, -.5, "Reward Probability of Stimuli", xpd=T, cex=1)

# Avoid Loss Context

plot_subset <- data_q[,c(which(names(data_q)=="Q_90_loss_corr"),which(names(data_q)=="Q_70_loss_corr"),which(names(data_q)=="Q_50_loss_corr"))]
expectationplot(plot_subset, ylab="action value estimates")

#add also the values of the incorrect choices
plot_subset <- data_q[,c(which(names(data_q)=="Q_90_loss_incorr"),which(names(data_q)=="Q_70_loss_incorr"),which(names(data_q)=="Q_50_loss_incorr"))]

for (c in 1:ncol(plot_subset)){
  lines(xaxis, plot_subset[,c], col=viridis_pal()(ncol(plot_subset))[c], lwd=2, lty=2) 
}

# insert dashed lines for programmed reward-probability
lines(x=c(0,480), y=c(0.5,0.5), col=viridis_pal()(3)[3], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.7,0.7), col=viridis_pal()(3)[2], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.9,0.9), col=viridis_pal()(3)[1], lty=2, lwd=1)


# Annotations
text(240, 1.1, "Avoid Loss", xpd=T, cex=1.2)

# Legend
legend(-100, -.5, c("90 %", "70 %", "50 %"),
       col=viridis_pal()(3), pch=19, bty="n", xpd=T,
       x.intersp=.7, ncol=5, lty=c(1)) # with default ncol=1, legend is vertical
dev.off()

# Solid lines indicate the estimates for the correct action in response to each respective stimulus while dashed indicate the incorrect one. Upper panels show stimuli learned in the reward context while lower panels show stimuli learned in the avoid punishment context.

#### Create plots 

# Separate plots for reward and loss context; separate lines for actions

# open a plot device
png("plots/action_values_full_range.png", width=7, height=3, unit="in", res=2400) #large enough

par(mfcol=c(1,2), mai=c(.7,.7,.3,.05), mgp=c(1.3,.4,0), tcl=-.25)

# First Plot: Get Reward Context

# subset data
plot_subset <- data_q[,c(which(names(data_q)=="Q_90_reward_corr"),which(names(data_q)=="Q_70_reward_corr"),which(names(data_q)=="Q_50_reward_corr"),which(names(data_q)=="Q_50_reward_incorr"),which(names(data_q)=="Q_70_reward_incorr"),which(names(data_q)=="Q_90_reward_incorr"))]

# create plot
expectationplot(plot_subset, ylab="action value estimates")

# insert dashed lines for programmed reward-probability
lines(x=c(0,480), y=c(0.1,0.1), col=viridis_pal()(6)[6], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.3,0.3), col=viridis_pal()(6)[5], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.5,0.5), col=viridis_pal()(6)[4], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.5,0.5), col=viridis_pal()(6)[3], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.7,0.7), col=viridis_pal()(6)[2], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.9,0.9), col=viridis_pal()(6)[1], lty=2, lwd=1)

# title above plot
text(240, 1.1, "Get Reward", xpd=T, cex=1.2)

# Legend (pt 1) below plot
text(60, -.5, "Reward/No Loss Probability", xpd=T, cex=1)
legend(210, -.58, c("90 %", "70 %"),
       col=viridis_pal()(6)[1:2], bty="n", xpd=T, # evtl fill = T mal probieren ## rausgenommen:  pch=19, 
       x.intersp=.7, ncol=5, lty=c(1), lwd=1.5) # with default ncol=1, legend is vertical

# Second plot: Avoid Loss Context

# subset data
plot_subset <- data_q[,c(which(names(data_q)=="Q_90_loss_corr"),which(names(data_q)=="Q_70_loss_corr"),which(names(data_q)=="Q_50_loss_corr"),which(names(data_q)=="Q_50_loss_incorr"),which(names(data_q)=="Q_70_loss_incorr"),which(names(data_q)=="Q_90_loss_incorr"))]

# create plot
expectationplot(plot_subset, ylab="action value estimates")

# insert dashed lines for programmed reward-probability
lines(x=c(0,480), y=c(0.1,0.1), col=viridis_pal()(6)[6], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.3,0.3), col=viridis_pal()(6)[5], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.5,0.5), col=viridis_pal()(6)[4], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.5,0.5), col=viridis_pal()(6)[3], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.7,0.7), col=viridis_pal()(6)[2], lty=2, lwd=1)
lines(x=c(0,480), y=c(0.9,0.9), col=viridis_pal()(6)[1], lty=2, lwd=1)

# Legend pt. 2
legend(-150, -.58, c("50 %", "50 %", "30 %", "10 %"),
       col=viridis_pal()(6)[3:6], bty="n", xpd=T, # evtl fill = T mal probieren ## rausgenommen:  pch=19, 
       x.intersp=.7, ncol=5, lty=c(1), lwd=1.5) # with default ncol=1, legend is vertical

# title above plot
text(240, 1.1, "Avoid Loss", xpd=T, cex=1.2)

dev.off()