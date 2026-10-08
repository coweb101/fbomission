### FBOmission - Characterization of apparent P3 waveform in omission trials

# (CW, 06/2026)

#### Routine & Info ####

## set working directory and read data
rm(list=ls())
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss")

#### Omission Data ####

data <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_700pre_1500post_anonym.csv") # read data again
data <- data.table::setDT(data) # convert to data.table
data <- subset(data, omission==1) # subset to trials with omitted feedback

segment_length <- 550 # samplepoints (-700 to 1500 ms with sampling rate 500 Hz)


line_variable <- data$valence # separate lines according to which variable?
# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data and get latencies and relative amplitudes for maxima for the frontocentral cluster ####


electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

data_avg_all <- data.frame()

# First step: Average across trials of each electrode separately for each participant and each 'unique condition' 

for (e in 1:length(electrodes)){
  
  electrode <- electrodes[e]
  
  print(electrode)
  
  # Subset data of electrode e
  data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
  
  # Ensure that amplitude data is numeric
  data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
  
  # Create averages for each participant separately for the unique condition combinations
  data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
  
  data_avg$electrode <- rep(electrode, times=nrow(data_avg))
  names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
  
  
  data_avg_all <- rbind(data_avg_all, data_avg)
  #data_se_all <- rbind(data_se_all, data_se)
}

data_avg_all <- do.call(cbind.data.frame, data_avg_all)

# Second step: Average across electrodes, and participants
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")] # Delete electrode column
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")] # Delete id column
data_avg <- c() #ensure this object is new
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean) # average across participants (stepwise procedure to get se for average across participants)

## Compute standard errors 

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)

## data_se and data_avg each have two rows (-1 and 1 for negative and positive feedback, respectively)
## data_se and data_avg each have 551 columns, whereas the first column indicates feedback valence (-1 and 1)

# transpose data
data_avg <- as.data.frame(t(data_avg))
data_se <- as.data.frame(t(data_se))

# rename columns and delete first row
names(data_se) <- c("negative", "positive")
data_se <- data_se[-1,]

names(data_avg) <- c("negative", "positive")
data_avg <- data_avg[-1,]

## latency
which(data_avg$positive == max(data_avg$positive))*4-704 
which(data_avg$negative == max(data_avg$negative))*4-704
# for both pos and neg, the most positive amplitude is observed at 444 ms

## amplitude

### based on the topographies and the ERPs, the deflection starts around 200 ms
### therefore, I take the time window between -200 to 200 as a baseline time
## window to extract a maximal amplitude relative to this reference

# -200 is samplepoint: 126 (because (-200+704)/4)
# 200 is samplepoint: 226

max(data_avg$positive) - mean(data_avg[126:226, "positive"])
max(data_avg$negative) - mean(data_avg[126:226, "negative"])

data_se[which(data_avg$positive == max(data_avg$positive)), "positive"]
data_se[which(data_avg$negative == max(data_avg$negative)), "negative"]

#### Aggregate data and get latencies and relative amplitudes for maxima for the centroparietal cluster ####

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")

data_avg_all <- data.frame()

# First step: Average across trials of each electrode separately for each participant and each 'unique condition' 

for (e in 1:length(electrodes)){
  
  electrode <- electrodes[e]
  
  print(electrode)
  
  # Subset data of electrode e
  data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
  
  # Ensure that amplitude data is numeric
  data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
  
  # Create averages for each participant separately for the unique condition combinations
  data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
  
  data_avg$electrode <- rep(electrode, times=nrow(data_avg))
  names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
  
  
  data_avg_all <- rbind(data_avg_all, data_avg)
  #data_se_all <- rbind(data_se_all, data_se)
}

data_avg_all <- do.call(cbind.data.frame, data_avg_all)

# Second step: Average across electrodes, and participants
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")] # Delete electrode column
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")] # Delete id column
data_avg <- c() #ensure this object is new
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean) # average across participants (stepwise procedure to get se for average across participants)

## Compute standard errors 

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)

## data_se and data_avg each have two rows (-1 and 1 for negative and positive feedback, respectively)
## data_se and data_avg each have 551 columns, whereas the first column indicates feedback valence (-1 and 1)

# transpose data
data_avg <- as.data.frame(t(data_avg))
data_se <- as.data.frame(t(data_se))

# rename columns and delete first row
names(data_se) <- c("negative", "positive")
data_se <- data_se[-1,]

names(data_avg) <- c("negative", "positive")
data_avg <- data_avg[-1,]

## latency (because values in the first 200ms are even more positive, I restrict this maximum search to timepoints after 0)
which(data_avg$positive == max(data_avg[175:550,"positive"]))*4-704 # most pos amp for pos fb at 528 ms
which(data_avg$negative == max(data_avg[175:550,"negative"]))*4-704 # most pos amp for neg fb at 664 ms

## amplitude

### based on the topographies and the ERPs, the deflection starts around 200 ms
### therefore, I take the time window between -200 to 200 as a baseline time
## window to extract a maximal amplitude relative to this reference

# -200 is samplepoint: 126 (because (-200+704)/4)
# 200 is samplepoint: 226

max(data_avg[175:550,"positive"]) - mean(data_avg[126:226, "positive"])
max(data_avg[175:550,"negative"]) - mean(data_avg[126:226, "negative"])

data_se[which(data_avg$positive == max(data_avg[175:550,"positive"])), "positive"]
data_se[which(data_avg$negative == max(data_avg[175:550,"negative"])), "negative"]

#### Display Data ####

data <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_700pre_1500post_anonym.csv") # read data again
data <- data.table::setDT(data) # convert to data.table
data <- subset(data, omission==-1) # subset to trials with displayed feedback

segment_length <- 550 # samplepoints (-700 to 1500 ms with sampling rate 500 Hz)


line_variable <- data$valence # separate lines according to which variable?
# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data and get latencies and relative amplitudes for maxima for the frontocentral cluster ####


electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

data_avg_all <- data.frame()

# First step: Average across trials of each electrode separately for each participant and each 'unique condition' 

for (e in 1:length(electrodes)){
  
  electrode <- electrodes[e]
  
  print(electrode)
  
  # Subset data of electrode e
  data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
  
  # Ensure that amplitude data is numeric
  data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
  
  # Create averages for each participant separately for the unique condition combinations
  data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
  
  data_avg$electrode <- rep(electrode, times=nrow(data_avg))
  names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
  
  
  data_avg_all <- rbind(data_avg_all, data_avg)
  #data_se_all <- rbind(data_se_all, data_se)
}

data_avg_all <- do.call(cbind.data.frame, data_avg_all)

# Second step: Average across electrodes, and participants
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")] # Delete electrode column
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")] # Delete id column
data_avg <- c() #ensure this object is new
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean) # average across participants (stepwise procedure to get se for average across participants)

## Compute standard errors 

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)

## data_se and data_avg each have two rows (-1 and 1 for negative and positive feedback, respectively)
## data_se and data_avg each have 551 columns, whereas the first column indicates feedback valence (-1 and 1)

# transpose data
data_avg <- as.data.frame(t(data_avg))
data_se <- as.data.frame(t(data_se))

# rename columns and delete first row
names(data_se) <- c("negative", "positive")
data_se <- data_se[-1,]

names(data_avg) <- c("negative", "positive")
data_avg <- data_avg[-1,]

## latency
which(data_avg$positive == max(data_avg$positive))*4-704 
which(data_avg$negative == max(data_avg$negative))*4-704
# for both pos and neg, the most positive amplitude is observed at 444 ms

## amplitude

### based on the topographies and the ERPs, the deflection starts around 200 ms
### therefore, I take the time window between -200 to 200 as a baseline time
## window to extract a maximal amplitude relative to this reference

# -200 is samplepoint: 126 (because (-200+704)/4)
# 200 is samplepoint: 226

max(data_avg$positive) - mean(data_avg[126:226, "positive"])
max(data_avg$negative) - mean(data_avg[126:226, "negative"])

data_se[which(data_avg$positive == max(data_avg$positive)), "positive"]
data_se[which(data_avg$negative == max(data_avg$negative)), "negative"]

#### Aggregate data and get latencies and relative amplitudes for maxima for the centroparietal cluster ####

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")

data_avg_all <- data.frame()

# First step: Average across trials of each electrode separately for each participant and each 'unique condition' 

for (e in 1:length(electrodes)){
  
  electrode <- electrodes[e]
  
  print(electrode)
  
  # Subset data of electrode e
  data_avg <- data[,c(which(colnames(data)=="id"), which(colnames(data)=="uniquecond"), which(substr(colnames(data),1,nchar(electrode))==electrode)), with=F]
  
  # Ensure that amplitude data is numeric
  data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode)] <- lapply(data_avg[,which(substr(colnames(data_avg),1,nchar(electrode))==electrode),with=F], as.numeric)
  
  # Create averages for each participant separately for the unique condition combinations
  data_avg <- aggregate(. ~ id + uniquecond, data= data_avg, mean)
  
  data_avg$electrode <- rep(electrode, times=nrow(data_avg))
  names(data_avg) <- c("id", "uniquecond", 1:segment_length, "electrode")
  
  
  data_avg_all <- rbind(data_avg_all, data_avg)
  #data_se_all <- rbind(data_se_all, data_se)
}

data_avg_all <- do.call(cbind.data.frame, data_avg_all)

# Second step: Average across electrodes, and participants
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")] # Delete electrode column
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")] # Delete id column
data_avg <- c() #ensure this object is new
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean) # average across participants (stepwise procedure to get se for average across participants)

## Compute standard errors 

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)

## data_se and data_avg each have two rows (-1 and 1 for negative and positive feedback, respectively)
## data_se and data_avg each have 551 columns, whereas the first column indicates feedback valence (-1 and 1)

# transpose data
data_avg <- as.data.frame(t(data_avg))
data_se <- as.data.frame(t(data_se))

# rename columns and delete first row
names(data_se) <- c("negative", "positive")
data_se <- data_se[-1,]

names(data_avg) <- c("negative", "positive")
data_avg <- data_avg[-1,]

## latency (because values in the first 200ms are even more positive, I restrict this maximum search to timepoints after 0)
which(data_avg$positive == max(data_avg[175:550,"positive"]))*4-704 # most pos amp for pos fb at 528 ms
which(data_avg$negative == max(data_avg[175:550,"negative"]))*4-704 # most pos amp for neg fb at 664 ms

## amplitude

### based on the topographies and the ERPs, the deflection starts around 200 ms
### therefore, I take the time window between -200 to 200 as a baseline time
## window to extract a maximal amplitude relative to this reference

# -200 is samplepoint: 126 (because (-200+704)/4)
# 200 is samplepoint: 226

max(data_avg[175:550,"positive"]) - mean(data_avg[126:226, "positive"])
max(data_avg[175:550,"negative"]) - mean(data_avg[126:226, "negative"])

data_se[which(data_avg$positive == max(data_avg[175:550,"positive"])), "positive"]
data_se[which(data_avg$negative == max(data_avg[175:550,"negative"])), "negative"]
