## FBOmission - Behavioural Data Analysis

## CW, last modified 10/2026

#### General settings & loading of packages ####
remove(list = ls()) # clear workspace
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss") # set wd if necessary

# for data wrangling
library(plyr)
library(reshape2)

# for statistical tests
library(lme4)
library(lmerTest)
library(emmeans)

# for plots
library(viridis)
library(ggeffects)

#general settings
options(width=500, scipen=6, digits=8)
options("jtools-digits" = 3)
control=glmerControl(optCtrl=list(maxfun=1e6))

sink("behavioural_analysis_output.txt")

set.seed(17)

#### Learning Performance (Choice Accuracy) ####

## Aggregation of behavioural data 

# In the analysis examining accuracy of choices (i.e. evidence of learning),
# only responses for the 70% and 90% stimuli are included as there was no right
# or wrong response for the 50% stimulus. Non-learners are defined as having 
# responded less accurate than chance and will be excluded from all further analyses.

# Read behavioural data

behav <- data.table::fread("aggregated_data/FBOmiss_behaviour_immediate_anonymized.csv")


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

# Invalid key presses (choice = 9): four of these trials were nevertheless
# followed by feedback (according to the log file); two of them were logged
# as incorrect (correct = 0). As they are no valid choices, "correct" is set to NA for all 
# of them (as in the RL models, see dataprep_PE_paper.m). This must be done
# before choice is recoded to NA.
behav[which(behav$choice == 9), "correct"] <- NA

# Recode "9" to NAs in "correct" and "choice" columns
behav[which(behav$correct == 9), "correct"] <- NA
behav[which(behav$choice == 9), "choice"] <- NA

# Four responses were logged marginally after the response deadline
# (2506-2517 ms) but were accepted by Presentation and followed by feedback;
# they are treated as valid responses (as in the RL models). Their logged
# response code is "incorrect" regardless of the key pressed, so accuracy is
# recoded from the chosen key and the participant's correct key for the
# respective stimulus (taken from trials with valid response codes).
# (This changes 1 of the 3 affected trials of stimuli 1-4 to correct.)

correct_key <- unique(behav[correct == 1 & response_time <= 2500 & stim %in% 1:4,
                            .(id, stim, correct_key = choice)])
stopifnot(!anyDuplicated(correct_key[, .(id, stim)])) # one correct key per stimulus
late <- which(behav$response_time > 2500 & !is.na(behav$choice) &
                behav$feedback %in% c(11, 1, -11, -1) & behav$stim %in% 1:4)
behav$correct[late] <- as.integer(
  behav$choice[late] == correct_key$correct_key[match(paste(behav$id[late], behav$stim[late]),
                                                      paste(correct_key$id, correct_key$stim))])

# Stimulus 5 & 6 had a reward probability of 50%; accordingly, there is no "correct choice".
# Therefore I recode the correct column for these stimuli before evaluating accuracy rates.
behav[which(behav$stim == 5|behav$stim == 6), "correct"] <- NA

# Note: Trials without a response within the deadline received no feedback
# (feedback = 9, context = NA) and are therefore not included in the GLMM.
# With the recoding above, the GLMM comprises 15,300 trials.

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
  behav <- behav[!(id %in% non_learner)]
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

# convert behav to data.frame which is needed for the indexing as done in the following loop
behav <- as.data.frame(behav)


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

behav_obs <- behav


#### Best BIC Model (2.1) ####
#### Analysis of accuracy data from simulated_data ####

behav_sim <- read.csv("aggregated_data/FBOmiss_behav_recovered_BIC.csv", sep=",", header=T)

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
png("plots/accuracy_posterior_prediction_BIC.png", width = 3.8, height = 1.9, units = "in", res = 2400)

par(mfrow = c(1, 2), mai = c(0.38, 0.4, 0.2, 0.05), mgp = c(1.4, 0.35, 0), tcl = -0.25, cex = 0.8)

for (ctx in c(1, -1)) { # loop through contexts
  
  # subset predicted and observed data
  plot_data <- plotdata[plotdata$context == ctx, ]
  obs_data <- agg_data[agg_data$context == ctx, ]
  
  ctx_label <- ifelse(ctx == 1, "Get Reward", "Avoid Loss")
  
  # create empty plot (for first plot with y-axis label)
  if (ctx==1){
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "Choice Accuracy (%)",
         xaxt = "n", yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
  } else{
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "", xaxt = "n",
         yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
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

#### Analysis of (Computationally-modelled) Learning Rates | BIC Model (M2.1) ####

## Data preparation

# read data
parameter_BIC <- read.csv("aggregated_data/FBOmiss_learning_parameter_BIC.csv")

# show descriptive statistics for parameters of rl model
params_BIC <- names(parameter_BIC)[-which(names(parameter_BIC)=="filename"|names(parameter_BIC)=="BIC")]
psych::describeBy(parameter_BIC[,params_BIC])

## reduce data (to learning rates, and filenames):
parameter_BIC <- parameter_BIC[,c(grep("filename", names(parameter_BIC)), grep("alpha", names(parameter_BIC)))]

# and reshape
data_learning_rates_BIC <- reshape2::melt(parameter_BIC, id = c("filename"))

# recode variables
data_learning_rates_BIC$learning_rate <- as.numeric(data_learning_rates_BIC$value)

data_learning_rates_BIC[grep("pos", data_learning_rates_BIC$variable), "valence"] <- "pos"
data_learning_rates_BIC[grep("neg", data_learning_rates_BIC$variable), "valence"] <- "neg"

## Transform data for approximating normal distribution

# cube-root transformation as for model 5b (chosen based on the normal qq-plots)
data_learning_rates_BIC$learning_original <- data_learning_rates_BIC$learning_rate
data_learning_rates_BIC$learning_rate <- data_learning_rates_BIC$learning_rate^(1/3)

# Plot distribution
qqnorm(data_learning_rates_BIC$learning_rate, pch = 1, frame = FALSE)
qqline(data_learning_rates_BIC$learning_rate, col = "steelblue", lwd = 2)

## Statistical analysis with linear mixed effects models: Data Preparation

data_BIC <- data_learning_rates_BIC

# code fixed effect of valence
data_BIC$valence <- ifelse(data_BIC$valence == "pos", 1, ifelse(data_BIC$valence == "neg", -1, NA))


## LME

# with two learning rates per participant, this model is equivalent to a paired t-test
model_learning_rates_BIC <- lmer(learning_rate ~ valence +
                                   (1|filename),
                                 data=data_BIC, REML=T,
                                 control=lmerControl(optimizer='bobyqa', optCtr = list(maxfun = 1e9)))
summary(model_learning_rates_BIC)
performance::r2(model_learning_rates_BIC,  tolerance = 1e-10)

# estimated marginal means (cube-root scale) and contrast positive vs. negative
emmeans::emmeans(model_learning_rates_BIC,
                 specs=pairwise~valence,
                 infer=TRUE,
                 lmer.df="asymp")

# check: paired t-test on the same (cube-root transformed) learning rates
t.test(parameter_BIC$alpha_pos^(1/3), parameter_BIC$alpha_neg^(1/3), paired = TRUE)

# inspect distribution of residuals
qqnorm(residuals(model_learning_rates_BIC))
qqline(residuals(model_learning_rates_BIC))

### Visualize

plotdata <- ggeffect(model_learning_rates_BIC, terms = c("valence"))

# Rename column and factorize
names(plotdata)[1] <- "valence"
plotdata$valence <- factor(plotdata$valence, levels = c(1, -1), labels = c("Positive", "Negative"))

# define colors (as for model 5b)
valence_colors <- c("blue", "green3")
names(valence_colors) <- c("Positive", "Negative")

# raw data
raw_subset <- data_BIC
raw_subset$valence <- factor(raw_subset$valence, levels = c(1, -1), labels = c("Positive", "Negative"))
raw_subset$x <- as.numeric(raw_subset$valence) +
  runif(nrow(raw_subset), -0.3, 0.3)  # jitter

# create plot for learning rates (arrows and bold points for transformed data as analyzed;
# and untransformed raw data points underneath for an impression of raw learning rate size)

png("plots/learning_rates_predicted_BIC.png", width = 1.9, height = 1.9, units = "in", res = 2400)

par(mfrow = c(1, 1), mai = c(0.38, 0.4, 0.3, 0.05), mgp = c(1.6, 0.35, 0), tcl = -0.25, cex = 0.8)

# Create empty plot
plot(1, type = "n", ylim = c(0, 1), xlim = c(0.5, 2.5),
     xaxt = "n", xlab = "", ylab = "Learning Rate",
     main = "", las = 1, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")

# Raw points
points(raw_subset$x, raw_subset$learning_original,
       pch = 16,
       col = adjustcolor(valence_colors[as.character(raw_subset$valence)], alpha.f = 0.4),
       cex = 0.5)

# Error bars
x_vals <- as.numeric(plotdata$valence)

arrows(x0 = x_vals, y0 = plotdata$conf.low, x1 = x_vals, y1 = plotdata$conf.high,
       angle = 90, code = 3, length = 0.05,
       col = valence_colors[as.character(plotdata$valence)])

# Predicted points
points(x_vals, plotdata$predicted,
       pch = 16,
       col = valence_colors[as.character(plotdata$valence)],
       cex = 0.8)

# Axis & label
axis(1, at = 1:2, labels = levels(plotdata$valence))
text(1.5, -0.31, "Feedback Valence", xpd = TRUE, cex = .9)

dev.off()

#### Comparison of fitted with recovered parameters ####

# load fitted (anonymized) parameter
parameter_fitted <- read.csv("aggregated_data/FBOmiss_learning_parameter_BIC.csv")

# load recovered parameter
parameter_recovered <- read.csv("aggregated_data/FBOmiss_learning_parameter_recovered_BIC.csv")

# create list with free parameter
list_parameter <- names(parameter_fitted)
list_parameter <- list_parameter[-c(which(list_parameter=="filename"), which(list_parameter=="BIC"))]

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
    
    if (list_parameter[j] != "beta") { 
      recovered_mat[, i] <- parameter_recovered[, colname]^(1/3) # cube-root transfrom learning rates
    } else {
      recovered_mat[, i] <- log(parameter_recovered[, colname]) # log-transform inverse temperature (right-skewed)
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

    if (!exists("plot_title")){
      plot_title <- paste0("Learning Rate (", valence, " Feedback)")
    }
    
  } else {
    fitted <- log(parameter_fitted[, list_parameter[j]]) # log-transform inverse temperature (right-skewed)
    xlim_use <- log(c(0.01, 100))
    ylim_use <- log(c(0.01, 100))
    if (!exists("plot_title")){
      plot_title <- "Inverse Temperature"
    }
  }
  
  # store correlation to write in plot title
  corr_coefficient <- cor(recovered_avg, fitted)
  correlations[1, list_parameter[j]] <- corr_coefficient
  
  
  # create plot
  png_filename <- file.path(paste0("plots/scatter_", list_parameter[j], "_BIC.png"))  
  png(png_filename, width = 1.9, height = 1.9, units = "in", res = 2400)
  
  par(mfrow = c(1, 1), mai = c(0.3, 0.35, 0.3, 0.25), mgp = c(1.6, 0.3, 0), tcl = -0.15, cex = 0.8)
  
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
    cex.axis = .7, las=1 , xaxt="n",
    yaxt = "n"
  )
  
  
  title(ylab="Fitted", line=1.3)
  title(xlab="Recovered", line=1.0)
  
  if (list_parameter[j] != "beta") {
    
    lr_ticks <- c(0, 0.01, 0.1, 0.5, 1)
    axis(1, at=lr_ticks^(1/3), labels=lr_ticks, las=1, cex.axis=.7, mgp=c(1.6,.001,0))
    axis(2, at=lr_ticks^(1/3), labels=lr_ticks, las=1, cex.axis=.7)
    
  }else{ # adjust axes for beta (log scale, labels on original scale)
    
    beta_ticks <- c(0.01, 0.1, 1, 10, 100)
    axis(1, at=log(beta_ticks), labels=beta_ticks, las=1, cex.axis=.7, mgp=c(1.6,.001,0))
    axis(2, at=log(beta_ticks), labels=beta_ticks, las=1, cex.axis=.7)
    
  }
  
  
  # add regression line
  abline(lm(fitted ~ recovered_avg), col = "darkgreen", lwd = 2)
  # add line for perfect correlation
  abline(a = 0, b = 1, col = "darkgray", lwd = 2, lty = 2)
  
  dev.off()
  
  print(plot_title)
  rm(plot_title)
  
}


#### Merge PE and behavioural data and visualize Action Value Estimates ####

# Read action values derived from reinforcement learning models based on the empirical data
pe_data <- read.csv("aggregated_data/FBOmiss_immediate_Q_values_and_PEs_BIC.csv")
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
data.table::fwrite(data,"aggregated_data/behav_and_pe_data_concatenated_BIC.csv", row.names=F) # save data in folder
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

# Create empty data frame with 12 columns (Q-values) and 480 rows (trials)
data_q <- data.frame(matrix(0, ncol = length(columns), nrow = ntrials)) # ncol = number of different Q-values, nrow = number of trials

# Loop through the Q-value columns and create averages of reward expectations
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
png("plots/action_values_correct_BIC.png", width=7, height=3, unit="in", res=2400) #large enough

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

# Solid lines indicate the estimates for the correct action in response to each respective stimulus while dashed indicate the incorrect one. Left panels show stimuli learned in the reward context while right panels show stimuli learned in the avoid punishment context.

#### Create plots 

# Separate plots for reward and loss context; separate lines for actions

# open a plot device
png("plots/action_values_full_range_BIC.png", width=7, height=3, unit="in", res=2400) #large enough

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


#### 2nd best BIC Model (5.1) ####
#### Analysis of accuracy data from simulated_data ####

behav_sim <- read.csv("aggregated_data/FBOmiss_behav_recovered_5b.csv", sep=",", header=T)

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

#### Resolve interactions

# Resolve interaction
emmeans::emtrends(accuracy_glmer_recovered, 
                  specs=pairwise~rew_prob| trialno_centered, 
                  var="trialno_centered",
                  infer=TRUE,
                  adj="bonf", lmer.df="asymp")
# Resolve interaction
emmeans::emtrends(accuracy_glmer_recovered, 
                  specs=pairwise~context| trialno_centered, 
                  var="trialno_centered",
                  infer=TRUE,
                  adj="bonf", lmer.df="asymp")
# Resolve interaction
emmeans::emmeans(accuracy_glmer_recovered, 
                 specs=pairwise~rew_prob| context , 
                 infer=TRUE,
                 adj="bonf", lmer.df="asymp")

emm2 <- emmeans::emmeans(accuracy_glmer_recovered, ~ rew_prob * context)
summary(contrast(emm2, interaction = c("pairwise", "pairwise")), infer = TRUE)

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
png("plots/accuracy_posterior_prediction_5b.png", width = 3.8, height = 1.9, units = "in", res = 2400)

par(mfrow = c(1, 2), mai = c(0.38, 0.4, 0.2, 0.05), mgp = c(1.4, 0.35, 0), tcl = -0.25, cex = 0.8)

for (ctx in c(1, -1)) { # loop through contexts
  
  # subset predicted and observed data
  plot_data <- plotdata[plotdata$context == ctx, ]
  obs_data <- agg_data[agg_data$context == ctx, ]
  
  ctx_label <- ifelse(ctx == 1, "Get Reward", "Avoid Loss")
  
  # create empty plot (for first plot with y-axis label)
  if (ctx==1){
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "Choice Accuracy (%)",
         xaxt = "n", yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
  } else{
    plot(NA, xlim = c(-1, 1), ylim = c(0, 100), xlab = "", ylab = "", xaxt = "n",
         yaxt="n", main = ctx_label, cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.8, bty="n")
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
parameter <- read.csv("aggregated_data/FBOmiss_learning_parameter_5b.csv")

# show descriptive statistics for parameters of rl model
params <- names(parameter)[-which(names(parameter)=="filename"|names(parameter)=="BIC")]
psych::describeBy(parameter[,params])

## reduce data (to learning rates, and filenames):
parameter <- parameter[,c(grep("filename", names(parameter)), grep("alpha", names(parameter)))]

# and reshape
data_learning_rates <- reshape2::melt(parameter, id = c("filename"))

# recode variables
data_learning_rates$learning_rate <- as.numeric(data_learning_rates$value)

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

# code fixed effect of valence
data$valence <- ifelse(data$valence == "pos", 1, ifelse(data$valence == "neg", -1, NA))
# and effect of appearance/omission
data$omission <- ifelse(data$type == "omission", 1, ifelse(data$type == "presentation", -1, NA))


## LME

model_learning_rates <- lmer(learning_rate ~ valence*omission + 
                               (1|filename),
                             data=data, REML=T, 
                             control=lmerControl(optimizer='bobyqa', optCtr = list(maxfun = 1e9)))
summary(model_learning_rates)
performance::r2(model_learning_rates,  tolerance = 1e-10) 

# Plot interaction (with odds logit)
emmeans::emmip(model_learning_rates,valence ~ omission, CIs = TRUE, cov.reduce = range, lmer.df="asymp")

# Resolve interaction
emmeans::emmeans(model_learning_rates, 
                 specs=pairwise~omission+valence, 
                 infer=TRUE,
                 adj="bonf", lmer.df="asymp")

# inspect distribution of residuals
qqnorm(residuals(model_learning_rates))
qqline(residuals(model_learning_rates))
# looks ok

### Visualize


plotdata <- ggeffect(model_learning_rates, terms = c("valence", "omission"))

# Rename column and factorize
names(plotdata)[1] <- "valence"
plotdata$valence <- factor(plotdata$valence, levels = c(1, -1), labels = c("Positive", "Negative"))

names(plotdata)[6] <- "omission"
plotdata$facet_group <- ifelse(plotdata$omission==-1,"Display","Omission")

# Split for facetting
facet_split <- split(plotdata, plotdata$facet_group)

# define colors
valence_colors <- c("blue", "green3")
names(valence_colors) <- c("Positive", "Negative")

# Plotting order
facet_order <- c("Display", "Omission")

# create plot for learning rates (arrows and bold points for transformed data as analyzed;
# and untransformed raw data points underneath for an impression of raw learning rate size)

png("plots/learning_rates_predicted_5b.png", width = 3.8, height = 1.9, units = "in", res = 2400)

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
raw_subset$valence <- factor(raw_subset$valence, levels = c(1, -1), labels = c("Positive", "Negative"))
raw_subset$omission <- ifelse(raw_subset$omission == -1, "Displayed", "Omitted")

raw_subset$x <- as.numeric(raw_subset$valence) +
  #ifelse(raw_subset$valence == "Positive", -0.15, 0.15) +
  runif(nrow(raw_subset), -0.3, 0.3)  # jitter

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
#df <- plotdata
x_vals <- as.numeric(df$valence)

arrows(x0 = x_vals, y0 = df$conf.low, x1 = x_vals, y1 = df$conf.high,
       angle = 90, code = 3, length = 0.05,
       col = valence_colors[as.character(df$valence)])

# Predicted points
points(x_vals, df$predicted,
       pch = 16,
       col = valence_colors[as.character(df$valence)],
       cex = 0.8)

# Axis & label
axis(1, at = 1:2, labels = levels(df$valence))
text(1.5, -0.31, "Feedback Valence", xpd = TRUE, cex = .9)

}

dev.off()


#### Comparison of fitted with recovered parameters ####

# load fitted (anonymized) parameter
parameter_fitted <- read.csv("aggregated_data/FBOmiss_learning_parameter_5b.csv")

# load recovered parameter
parameter_recovered <- read.csv("aggregated_data/FBOmiss_learning_parameter_recovered_5b.csv")

# create list with free parameter
list_parameter <- names(parameter_fitted)
list_parameter <- list_parameter[-c(which(list_parameter=="filename"), which(list_parameter=="BIC"))]

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
    
    if (list_parameter[j] != "beta") { 
      recovered_mat[, i] <- parameter_recovered[, colname]^(1/3) # cube-root transform learning rates
    } else {
      recovered_mat[, i] <- log(parameter_recovered[, colname]) # log-transform inverse temperature (right-skewed)
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
    appearance <- ifelse(grepl("omitted",list_parameter[j]), "Omission", "Display")
    
    if (!exists("plot_title")){
      plot_title <- paste0("LR (", valence, ", ", appearance, ")")
    }
    
  } else {
    fitted <- log(parameter_fitted[, list_parameter[j]]) # log-transform inverse temperature (right-skewed)
    xlim_use <- log(c(0.01, 100))
    ylim_use <- log(c(0.01, 100))
    if (!exists("plot_title")){
      plot_title <- "Inverse Temperature"
    }
  }
  
  # store correlation to write in plot title
  corr_coefficient <- cor(recovered_avg, fitted)
  correlations[1, list_parameter[j]] <- corr_coefficient
  
  # create plot
  png_filename <- file.path(paste0("plots/scatter_", list_parameter[j], "_5b.png"))  
  png(png_filename, width = 1.9, height = 1.9, units = "in", res = 2400)
  
  par(mfrow = c(1, 1), mai = c(0.3, 0.35, 0.3, 0.25), mgp = c(1.6, 0.3, 0), tcl = -0.15, cex = 0.8)
  
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
    cex.axis = .7, las=1 , xaxt="n",
    yaxt = "n" 
  )
  
  
  title(ylab="Fitted", line=1.3)
  title(xlab="Recovered", line=1.0)
  
  if (list_parameter[j] != "beta") {
    
    lr_ticks <- c(0, 0.01, 0.1, 0.5, 1)
    axis(1, at=lr_ticks^(1/3), labels=lr_ticks, las=1, cex.axis=.7, mgp=c(1.6,.001,0))
    axis(2, at=lr_ticks^(1/3), labels=lr_ticks, las=1, cex.axis=.7)
    
  }else{ # adjust axes for beta (log scale, labels on original scale)
    
    beta_ticks <- c(0.01, 0.1, 1, 10, 100)
    axis(1, at=log(beta_ticks), labels=beta_ticks, las=1, cex.axis=.7, mgp=c(1.6,.001,0))
    axis(2, at=log(beta_ticks), labels=beta_ticks, las=1, cex.axis=.7)
    
  }
  
  
  # add regression line
  abline(lm(fitted ~ recovered_avg), col = "darkgreen", lwd = 2)
  # add line for perfect correlation
  abline(a = 0, b = 1, col = "darkgray", lwd = 2, lty = 2)
  
  dev.off()
  
  print(plot_title)
  rm(plot_title)
  
}



#### Additional Analysis: Recovery of the difference between learning rates for displayed and omitted outcomes (M5.1) ####

# difference displayed - omitted (cube-root scale), for positive and negative outcomes
for (v in c("pos", "neg")) {
  
  # fitted difference per participant
  fitted_diff <- parameter_fitted[, paste0("alpha_", v, "presented")]^(1/3) -
    parameter_fitted[, paste0("alpha_", v, "omitted")]^(1/3)
  
  # recovered difference per participant: computed per simulation, then averaged across the 25 simulations
  recovered_diff_mat <- matrix(NA, nrow = nrow(parameter_recovered), ncol = 25)
  for (i in 1:25) {
    sim <- sprintf("%02d", i)
    recovered_diff_mat[, i] <- parameter_recovered[, paste0("alpha_", v, "presented_", sim)]^(1/3) -
      parameter_recovered[, paste0("alpha_", v, "omitted_", sim)]^(1/3)
  }
  recovered_diff <- rowMeans(recovered_diff_mat)
  
  # correlation and sign agreement
  print(paste0(v, " outcomes: r = ", round(cor(fitted_diff, recovered_diff), 2),
               "; same sign in ", sum(sign(fitted_diff) == sign(recovered_diff)),
               " of ", length(fitted_diff), " participants"))
  print(paste0("exceptions (fitted): ", paste(round(fitted_diff[sign(fitted_diff) != sign(recovered_diff)], 3), collapse = ", ")))
  print(paste0("mean fitted = ", round(mean(fitted_diff), 2), ", mean recovered = ", round(mean(recovered_diff), 2)))
}

#### Merge PE and behavioural data and visualize Action Value Estimates ####

# Read action values derived from reinforcement learning models based on the empirical data
pe_data <- read.csv("aggregated_data/FBOmiss_immediate_Q_values_and_PEs_5b.csv")
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
data.table::fwrite(data,"aggregated_data/behav_and_pe_data_concatenated_5b.csv", row.names=F) # save data in folder
pe_data <- data # rename in the environment for later merging of scripts with eeg_data

#### Expectation Plots 

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

# Create empty data frame with 12 columns (Q-values) and 480 rows (trials)
data_q <- data.frame(matrix(0, ncol = length(columns), nrow = ntrials)) # ncol = number of different Q-values, nrow = number of trials

# Loop through the Q value columns and create averages of reward expectations
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
png("plots/action_values_correct_5b.png", width=7, height=3, unit="in", res=2400) #large enough

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

# Solid lines indicate the estimates for the correct action in response to each respective stimulus while dashed indicate the incorrect one. Left panels show stimuli learned in the reward context while right panels show stimuli learned in the avoid punishment context.

#### Create plots 

# Separate plots for reward and loss context; separate lines for actions

# open a plot device
png("plots/action_values_full_range_5b.png", width=7, height=3, unit="in", res=2400) #large enough

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





#### Prepare table for modelfits ####


# read fits
fits <- data.table::fread("aggregated_data/modelfits.csv", header=T, sep = ",", dec = ".")

# to copy in word table
knitr::kable(
  as.data.frame(
    round(fits$AIC,2)) # AIC column
)

knitr::kable(
  as.data.frame(
    round(fits$BIC,2)) # BIC column
)

fits[which(fits$BIC==min(fits$BIC)),] # show model with lowest BIC
fits[order(fits$BIC),] # show all models ordered by BIC


#### Comparison of the two best models (M2.1 = model 2b, M5.1 = model 5b) ####

## Mean and summed differences (modelfits.csv contains means across participants)

# per-participant BIC values of both models
bic_2b <- read.csv("aggregated_data/FBOmiss_learning_parameter_BIC.csv")
bic_5b <- read.csv("aggregated_data/FBOmiss_learning_parameter_5b.csv")
stopifnot(all(bic_2b$filename == bic_5b$filename)) # same participant order in both files

n_participants <- nrow(bic_2b)

# BIC difference (M5.1 - M2.1): mean per participant and summed across participants
delta_bic_mean <- fits[model == "5b", BIC] - fits[model == "2b", BIC]
print(paste0("BIC (M5.1 - M2.1): mean = ", round(delta_bic_mean, 2),
             ", summed = ", round(delta_bic_mean * n_participants, 1)))

# AIC difference (M2.1 - M5.1): mean per participant and summed across participants
delta_aic_mean <- fits[model == "2b", AIC] - fits[model == "5b", AIC]
print(paste0("AIC (M2.1 - M5.1): mean = ", round(delta_aic_mean, 2),
             ", summed = ", round(delta_aic_mean * n_participants, 1)))


## Per-participant differences

# per-participant AIC from BIC: AIC = BIC - k*ln(n) + 2k
# (k = 3 for M2.1 and 5 for M5.1; n = trials with valid response and feedback, as in the model fitting)
behav_raw <- data.table::fread("aggregated_data/FBOmiss_behaviour_immediate_anonymized.csv")
n_valid <- behav_raw[choice %in% c(1, 2) & feedback %in% c(11, 1, -11, -1), .N, by = id]
n_valid <- n_valid$N[match(bic_2b$filename, n_valid$id)]
aic_2b <- bic_2b$BIC - 3 * log(n_valid) + 2 * 3
aic_5b <- bic_5b$BIC - 5 * log(n_valid) + 2 * 5

delta_bic <- bic_5b$BIC - bic_2b$BIC # positive = M2.1 better
delta_aic <- aic_2b - aic_5b         # positive = M5.1 better

# check: must equal the summed AIC difference above (363.7)
print(paste0("check summed AIC difference: ", round(sum(delta_aic), 1)))

# number of participants favoring each model
print(paste0("lower BIC for M2.1 in ", sum(delta_bic > 0), " of ", n_participants, " participants"))
print(paste0("lower AIC for M5.1 in ", sum(delta_aic > 0), " of ", n_participants, " participants"))

# distribution of BIC differences
print(paste0("delta BIC across all participants: median = ", round(median(delta_bic), 2),
             ", IQR = ", round(quantile(delta_bic, .25), 2), " to ", round(quantile(delta_bic, .75), 2)))


## Participants excluded from the EEG analyses

# same criterion as in the EEG scripts: SD of |PE| < .05 in display or omission trials
pe_data_2b <- data.table::fread("aggregated_data/behav_and_pe_data_concatenated_BIC.csv")
pe_sd <- pe_data_2b[!is.na(PE) & omission %in% c(-1, 1),
                    .(sd_abs_pe = sd(abs(PE))), by = .(id, omission)]
low_var_ids <- unique(pe_sd[sd_abs_pe < .05, id])
excluded <- bic_2b$filename %in% low_var_ids

# values of the excluded participants and their rank among all participants (1 = strongest preference for M2.1)
print(data.frame(id = bic_2b$filename[excluded],
                 delta_bic = round(delta_bic[excluded], 2),
                 delta_aic = round(delta_aic[excluded], 2),
                 rank_delta_bic = rank(-delta_bic)[excluded]))

# model comparison without the excluded participants
print(paste0("without excluded participants (n = ", sum(!excluded), "): ",
             "BIC mean = ", round(mean(delta_bic[!excluded]), 2),
             ", summed = ", round(sum(delta_bic[!excluded]), 1),
             "; lower BIC for M2.1 in ", sum(delta_bic[!excluded] > 0), " of ", sum(!excluded),
             "; lower AIC for M5.1 in ", sum(delta_aic[!excluded] > 0), " of ", sum(!excluded),
             "; AIC mean = ", round(mean(delta_aic[!excluded]), 2),
             ", summed = ", round(sum(delta_aic[!excluded]), 1)))


#### Stop saving output ####
sink()
