#### FBOMISS - EEG Data Aggregation

# CW,  last modified 10/2026

#### Info ####

# First: EEG data is read (BVA exportfiles)
# Second: Empty rows are added for excluded trials (BVA artifact rejection)
# or missing trials and the data set is merged with the havioural/pe-data.


#### Routine ####

remove(list = ls()) # clear workspace
getwd() # show current working directory
setwd("//psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss")


##### Read eeg and marker data (BVA export) ####

# where is the data (in relation to working directory)?
eegdatafolder <- "eeg_data/immediate/export_revision"
## when it is in a subfolder, e.g. "Raw Data EEG/Export"
## this folder should contain single-trial .dat and .vmrk files

# define the naming scheme of the export files (in this case, the filenames are composed of the ids and the string "_preprocessed")
naming_scheme <- "_preprocessed_no_baseline_extended_segment.dat"

segment_length <- 550 # no of sample points per segment

# List all files with specified pattern in folder "Export" and save them as list "datfiles"
datfiles <- list.files(eegdatafolder, pattern = naming_scheme)

#x <- datfiles[4] # to try commands inside of lapply with one file
#datfiles <- datfiles[-3]

# Apply the following to each entry of list "datfiles"
erp_data <- lapply(datfiles, function(x) {
  
  ### Collect basic data from current filename x (eeg-export filename) and filepath
  
  # subject code
  id <- substr(x, 1, 6) # first six characters of eeg-export-filename
  
  ### Read markerfile
  
  # In order to find respective markerfile:
  mrkname <- sub(".dat$",".vmrk",x) # replace ".dat" at the end ($) with .vmrk
  file <- paste0(eegdatafolder, "/", mrkname, collapse="") # paste full path+filename for .vmrk
  mrk <- readLines(file) # read markerfile
  
  # Find number of comment lines above marker info
  skip <- grep("Mk1", mrk)[1]-1 # grep() finds string in file (here:mrk) & returns number of first line with [1]
  mrk <- read.csv(file,skip=skip,header=FALSE) # use skip option to skip comment lines
  
  # Get Marker
  stim <- mrk[grep("Mk[0-9]+=Stimulus",mrk[,1]),] # stim <- mrk but only lines with Stimulus (in Index: find Mk + number + = Stimulus in erster Spalte vom mrk)
  stim <- sub("S","",stim[,2]) # substitutes "S" with empty string i second column
  stim <- as.numeric(stim)
  stim <- subset(stim, stim>234) # delete markers below 235 (which are no feedbackmarkers but may be in segments)
  
  ### Read the respective datfile & extract amplitudes of electrodes of interest
  ### and bring data for each trial in a respective row
  
  dat <- data.table::fread(paste0(eegdatafolder, "/", x, collapse=""), header=T, sep = " ", dec = ",") # replaced read.delim with fread//20240814
  dat <- as.data.frame(apply(dat, 2, as.numeric))# in case R classifies variables otherwise
  
  # create vector with names of all electrodes in the export file
  export_electrodes <- names(dat)
  
  # create empty dataframe
  df <- data.frame(matrix(nrow = length(stim), ncol =segment_length*length(export_electrodes))) 
  
  for (i in 1:length(export_electrodes)){
    
    current_electrode <- export_electrodes[i]
    
    amplitude_data <- dat[,current_electrode] # take data from electrode of interest
    amplitude_data <- as.data.frame(split(amplitude_data, 1:segment_length)) # split vector in trials (based on number of samplepoints per trial/segment_length)
    
    df[,((i-1)*segment_length+1):(i*segment_length)] <- amplitude_data # paste in large dataframe
    names(df)[((i-1)*segment_length+1):(i*segment_length)] <- paste0(current_electrode,"_", 1:segment_length) # name columns
  }
  
  id <- rep(id, times=length(stim))
  
  ### Return everything that should be saved in the resulting dataframe
  return(cbind.data.frame(id, stim, df))
  gc()
})

eeg_data <- do.call(rbind.data.frame, erp_data)

data.table::fwrite(eeg_data, "aggregated_data/bva_export_files_concatenated_700pre_1500post.csv", row.names=F) # save data

#### Number trials (for a later check on preserved order) ####

data <- eeg_data

trial <- c() # create empty trial vector

for (i in 1:length(unique(data$id))) {
  
  # build vector of 1 to number of rows the current filename has in data
  trial_temp <- c(1:nrow(data[which(data$id == unique(data$id)[i]),]))
  trial <- c(trial, trial_temp) # append to trial vector
  
}

data$trial <- trial # append trial vector to data

##### Insert empty rows for trials which were rejected during Artifact Scan in BVA #####

## Because we want to match the data with trial-by-trial PEs created in Matlab,
## we need to add empty rows for all trials that are not included in the EEG 
## data export

# List removed segments from artifact rejection (two files as there were
# separate preprocessing trees for experiment version a and b)
art <- readLines(paste0(eegdatafolder, "/", "Report_Artifact Rejection_Long.txt"))

seglist <- list()
j=1 # start counter
for (i in 1:length(grep("History File:", art))) { 
  filename <- art[grep("History File:", art)[i]+1]
  id <- substr(filename, 1,6) # first six characters of line
  
  startseg <- grep("The following segments have been removed:", art)[i]+1
  endseg <- grep("Artifact Type", art)[j]-6
  seg <- art[startseg:endseg] 
  seg <- as.numeric(unlist(strsplit(seg, split=", ")))
  
  seglist[[i]] <- cbind(id, seg)
  
  j=j+2
  # add 2 to counter (to grep the correct line in which the string "Artifact
  # Type" occurs the second time for the current subject; in contrast to the
  # strings "History File" and "... have been removed:", "Artifact Type" occurs 
  # 2 times in every subject report, therefore, the separate counter)
  
}

## Define function to add empty rows for removed segments
insertRows <- function (dataframe, newrows, index) {
  temp <- dataframe
  for (i in 1:nrow(newrows)){
    
    if (index[i]+i-1 == nrow(temp)+1 ) {
      temp <- rbind(temp, newrows[i,])
      # special case when index is last row which does not exist yet
    } else {
      temp <- as.data.frame(temp,stringsAsFactors=FALSE)
      # inserts empty row at position i of index
      temp[seq(index[i] + i,nrow(temp)+1),] <- temp[seq(index[i] + i - 1,nrow(temp)),]
      temp[index[i] + i - 1,] <- newrows[i,]
    }
  }
  
  return(temp)
}


# Add removed segments from Artifact Rejection
columnsno <- ncol(data)
insertdata <- data.frame() # New data frame

for (i in 1:length(seglist)) {
  x <- seglist[[i]]
  if (length(x) > 1) {
    
    gc() # garbage collection as the script does not run due to ram
    
    nnewrows <- nrow(x) # number new rows
    rows <- matrix(NA, ncol=columnsno, nrow=nnewrows)
    rows[,which(colnames(data)=="id")] <- x[,"id"]
    
    temp <- subset(data, id==unique(x[,"id"]))
    index <- as.numeric(x[,"seg"])
    index <- sort(index)
    index <- index - 0:(length(index)-1)
    tempdata <- insertRows(temp, rows, index)
    insertdata <- rbind(insertdata, tempdata)
    
    id <- x[1,"id"]
    print(paste0(nnewrows, " row(s) added to data of ", id))
  } else insertdata <- rbind(insertdata, subset(data, id==x[,"id"]))
}

table(insertdata$id) # check whether trials are now complete for each subject

data <- insertdata
remove(insertdata)

data.table::fwrite(data,"aggregated_data/bva_export_files_concatenated_filled_rows_700pre_1500post.csv", row.names=F) # save data
#data <- data.table::fread("aggregated_data/bva_export_files_concatenated_filled_rows_700pre_1500post.csv", header=T, sep = ",", dec = ".") 

##### Check whether intial order of trials was preserved #####

for (i in 1:length(unique(data$id))) {
  
  x <- data[which(data$id == unique(data$id)[i]), ]
  x <- as.data.frame(x)
  print(!is.unsorted(as.numeric(x[,"trial"]), na.rm=T))
  # print true if the order of trial is ascending
  # i.e. if only TRUE appears, everthing's in order
  
} # ok, cool


#### Anonymize EEG data ####

# get list created in Matlab at first point of anonymization (and which was also used in dataprep of accuracy data)
id_key <- data.table::fread("behavioural_data/immediate/interim_datasets/filenames_old_new.csv")

# replace pseudonymized code with anonymous id

for (i in 1:length(unique(id_key$original_id))){
  
  id_current <- unique(id_key$original_id)[i]
  data[which(data$id== id_current), "id"] <- id_key[which(id_key$original_id==id_current),2]
  
}

rm(id_key)
rm(id_current)


##### Read behav/PE-data #####

# read behavioural and PE data file created in
# FBOmission_Behavioural_Data_Analysis_anonym.R
pe_data <- data.table::fread("aggregated_data/behav_and_pe_data_concatenated.csv", header=T, sep = ",", dec = ".")


##### Check whether experimental order is still accurate in both data sets ####

# feedback coding in behavioural/pe data:
#   positive outcomes:
#     11 - feedback is presented which indicates a monetary gain
#     1 - feedback is omitted which indicates that a monetary loss was avoided
#   negative outcomes:
#     -1 - feedback is omitted which indicates that a monetary gain was missed
#     -11 - feedback is presented which indicates a monetary loss

# feedback coding in eeg data (first two numbers are either 23, 24 or 25 
# depending on the version of the experiment (a or b)), last number indicates 
# type and valence:
#   positive outcomes:
#     0 - feedback is presented which indicates a monetary gain
#     8 - feedback is omitted which indicates that a monetary loss was avoided
#   negative outcomes:
#     7 - feedback is omitted which indicates that a monetary gain was missed
#     9 - feedback is presented which indicates a monetary loss
#   other:
#     5 - no response

# code variable to compare in both data sets
pe_data$eeg_marker <- ifelse(pe_data$feedback==11, 0,
                             ifelse(pe_data$feedback==1, 8,
                                    ifelse(pe_data$feedback==-1, 7,
                                           ifelse(pe_data$feedback==-11, 9, NA))))
data$eeg_marker <- as.numeric(substr(data$stim,3,3))

# compare this column between data sets (row index created with order() sorts
# the column by id) and print the number of non-matching entries:
sum(pe_data[order(as.numeric(as.factor(pe_data$id))),"eeg_marker"] !=
      data[order(as.numeric(as.factor(data$id))),"eeg_marker"], na.rm=T)

##### Merge behavioural/pe- and eeg-data ##### 

# as we just checked and confirmed the accurate order of trials in both data
# sets via the marker columns, we can now add a new trial variable to both data
# sets which we the use to merge the data (to ensure that corresponding trials
# are matched for each participant):

pe_data$trial <- rep(1:480, length(unique(pe_data$id)))
data$trial <- rep(1:480, length(unique(data$id))) # initial variable "trial" (for order check, see above) is overwritten

library(plyr) # to use join() (merging data while preserving row order (merge() does that not...))

data <- plyr::join(pe_data, data, by = c("id", "trial"))


data.table::fwrite(data, file="aggregated_data/eeg_and_pe_data_merged_700pre_1500post_anonym.csv", row.names=F) # wanna save data in between?
#data <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_700pre_1500post_anonym.csv")


