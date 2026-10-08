### Feedback Omission - Trial Count Table

### First - Get number of trials with valid responses (within the response time window)

remove(list = ls()) # clear workspace
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss")
data <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_700pre_1500post_anonym.csv") # read concatenated data

data$valid_trial <- ifelse(!is.na(data$F7_1), 1, 0)

data <- as.data.frame(data)
trial_count <- aggregate(data$valid_trial, by=list(data$id, data$omission, data$valence), FUN="sum")

aggregate(trial_count$x, by=list(trial_count$Group.2, trial_count$Group.3), FUN="mean")
# Group.1 Group.2         x
# 1      -1      -1 100.08333
# 2       1      -1  99.58333
# 3      -1       1 134.91667
# 4       1       1 134.00000

aggregate(trial_count$x, by=list(trial_count$Group.2, trial_count$Group.3), FUN="range")
# Group.1 Group.2 x.1 x.2
# 1      -1      -1  71 132
# 2       1      -1  71 139
# 3      -1       1  99 162
# 4       1       1 104 164

## Group.1 --> Omission with -1 for Feedback Display, and 1 for Feedback Omission
## Group.2 --> Feedback Valence with -1 for Negative Feedback, and 1 for Positive Feedback
