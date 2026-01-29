#### FBOmission Multitemporal Results Plots ####
# CW (12/2025)

#### Routine & Info ####

## set working directory and read data
rm(list=ls())
getwd() # show current working directory
setwd("\\\\psychologie.ad.hhu.de/biopsych_experimente/Studien_Daten/2024_CB_CW_FBOmiss")

# required packages: scales, viridis, data.table, yarrr (functions are directly called when used via ::)


#### Omission | All Fixed Effect Coefficients | -700 to 1500 ms || Read data and set general settings ####

pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_omission_trials.csv", quote="") # read data again
pe_coefficients_centroparietal <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_omission_trials.csv", quote="") # read data again
data_grand_averages <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_7001500_anonym.csv", quote="") # read data again

segment_start <- 700
segment_end <- 1500

segment_length <- 550

# axis labels
axis_start <- 500
axis_label_distance <- 500

# axis ticks (in samplepoints)
axis_start_sample_points <- 50
axis_distance <- 125
zero_line_sample_points <- 175

ylim <- c(4,-4) # limits of y-axis in all plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/4)
x_lab_text <- 275


#### Open plot device ####
png(filename = paste0("plots/multitemp_results_omission_trials_all_fixed_effects_extended_segment_",as.character(Sys.Date()),".png"),
    width=8, height=7, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(3,2), mai=c(.43,.5,.4,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

#### Plot 1 (topleft): Effects at Frontocentral Cluster | PE ####

# get/rename data
pe_coefficients <- pe_coefficients_frontocentral
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction
significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot

## Plot significance indicators

# to loop through fixed effects:
effects <- c("significanteffect_pe", "significanteffect_valence", "significanteffect_interaction")

# Define the mixed color
mid_color <- scales::alpha("#006780", 0.7)
col <- c("purple", mid_color, "cornflowerblue")
pch <- c(15, 15, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}


## Add indicators for model convergence and singular fit

# add yellow circles for singular fit models
if (length(singular)>0){
  for (i in 1:length(singular)){
    points(singular[i], (ylim[1]-(.0375*abs(diff(ylim)))), pch=19, col="coral", xpd=T)
  }
}

# add red circles for non-converged model (overwrite singular fit if it is the same)
if (length(non_convergence)>0){
  for (i in 1:length(non_convergence)){
    points(non_convergence[i], ylim[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
  }
}
## Draw lines for coefficients

# main effect valence
lines(pe_coefficients$time, pe_coefficients$coef_valence*2, col=mid_color, lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_valence*2 + pe_coefficients$se_valence, 
                       rev(pe_coefficients$coef_valence*2 - pe_coefficients$se_valence))), col=yarrr::transparent(mid_color, trans.val = .7),xpd=T, border=NA)

# interaction valence + pe
lines(pe_coefficients$time, pe_coefficients$coef_interaction_pe_valence, col="cornflowerblue", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_interaction_pe_valence + pe_coefficients$se_interaction_pe_valence, 
                       rev(pe_coefficients$coef_interaction_pe_valence - pe_coefficients$se_interaction_pe_valence))),
        col=yarrr::transparent("cornflowerblue", trans.val = .7),xpd=T, border=NA)

# main effect pe
lines(pe_coefficients$time, pe_coefficients$coef_pe, col="purple", lwd = 3, xpd=T)
 
polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                       rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent("purple", trans.val = .7),xpd=T, border=NA)


## Annotations
# Axes

axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# y-axis label
title(ylab = expression(beta~"-coefficient"))

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset



#### Plot 2 (middleleft): Grand Averages Separately for Valence ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data

data <- subset(data, omission==1) # only trials with omitted feedback

# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)
# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$valence # separate lines according to which variable?

# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- subset(data_avg, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp


# standard error
temp <- subset(data_se, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


short_pos <- cbind.data.frame(short_pos, short_neg)

short_pos_se <- cbind.data.frame(short_pos_se, short_neg_se)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols=c("blue", "green3"), 
      lwd=3, 
      colsterror=c(yarrr::transparent("blue", trans.val = .7), yarrr::transparent("green3", trans.val = .7))) # topleft
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 3 (bottomleft): Grand Averages Separately for PEs ####

## get data
data <- aggr_data 
data <- data.table::setDT(data) # convert to data.table
data <- subset(data, omission==1) # subset to trials with omitted feedback

## define custom function

x <- seq(-segment_start, segment_end, length.out=segment_length) # define x-axis
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 
data$uniquecond <- data$pe_c_cat # rename to uniquecond to recyclue previous script

## choose electrodes and pipe general info
electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

## aggregate data 

data_avg_all <- data.frame()

# average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- data_avg
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- data_se
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)



gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols= viridis::viridis_pal()(ncol(short_pos)), 
      lwd=3, 
      colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 4 (topright): Effects at Centroparietal Cluster | PE ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot

## Plot significance indicators

# to loop through fixed effects:
effects <- c("significanteffect_pe", "significanteffect_valence", "significanteffect_interaction")

# Define the mixed color
mid_color <- scales::alpha("#006780", 0.7)
col <- c("purple", mid_color, "cornflowerblue")
pch <- c(15, 15, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}

## Add indicators for model convergence and singular fit

# add yellow circles for singular fit models
if (length(singular)>0){
  for (i in 1:length(singular)){
    points(singular[i], (ylim[1]-(.0375*abs(diff(ylim)))), pch=19, col="coral", xpd=T)
  }
}

# add red circles for non-converged model (overwrite singular fit if it is the same)
if (length(non_convergence)>0){
  for (i in 1:length(non_convergence)){
    points(non_convergence[i], ylim[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
  }
}


## Draw lines for coefficients

lines(pe_coefficients$time, pe_coefficients$coef_valence*2, col=mid_color, lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_valence*2 + pe_coefficients$se_valence, 
                       rev(pe_coefficients$coef_valence*2 - pe_coefficients$se_valence))), col=yarrr::transparent(mid_color, trans.val = .7),xpd=T, border=NA)


lines(pe_coefficients$time, pe_coefficients$coef_interaction_pe_valence, col="cornflowerblue", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_interaction_pe_valence + pe_coefficients$se_interaction_pe_valence, 
                       rev(pe_coefficients$coef_interaction_pe_valence - pe_coefficients$se_interaction_pe_valence))),
        col=yarrr::transparent("cornflowerblue", trans.val = .7),xpd=T, border=NA)


# 
# lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
# lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_pe, col="purple", lwd = 3, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
#                        rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
#                        rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)
# 
# 
polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                       rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent("purple", trans.val = .7),xpd=T, border=NA)


## Annotations
# Axes
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset


#### Plot 5 (middleleft): Grand Averages Separately for Valence ####

# get data
data <- aggr_data 
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data

#### Plot 2 & 3 (left middle and bottom): Grand Averages | Omitted | PE | Open plot & set global parameters 


# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)


# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


#### Prepare Grand Averages for Frontocentral Electrodes

#data <- aggr_data
data <- subset(data, omission==1) # only trials with omitted feedback

##### Categorize continuous pe variable for separate lines

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$valence # separate lines according to which variable?
#within_plot_variable <- data$valence # separate plots according to which variable?



# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- subset(data_avg, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp


# standard error
temp <- subset(data_se, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


short_pos <- cbind.data.frame(short_pos, short_neg)

short_pos_se <- cbind.data.frame(short_pos_se, short_neg_se)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols=c("blue", "green3"), 
      lwd=3, 
      colsterror=c(yarrr::transparent("blue", trans.val = .7), yarrr::transparent("green3", trans.val = .7))) # topleft
#title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 6 (bottomleft): Grand Averages Separately for PEs ####

## get data
data <- aggr_data
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data
data <- subset(data, omission==1) # subset to trials with omitted feedback

## define custom function

x <- seq(-segment_start, segment_end, length.out=segment_length) # define x-axis
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 
data$uniquecond <- data$pe_c_cat # rename to uniquecond to recyclue previous script

## choose electrodes and pipe general info
electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

## aggregate data 

data_avg_all <- data.frame()

# average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- data_avg
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- data_se
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)



gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols= viridis::viridis_pal()(ncol(short_pos)), 
      lwd=3, 
      colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
#title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Close plot ####
dev.off() # close device








#### Clear workspace ####
rm(list=ls())
#### Omission | All Fixed Effect Coefficients | -200 to 800 ms || Load data, set general settings and open plot ####

pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre200_post800_frontocentral_cluster_omission_trials.csv", quote="") # read data again
pe_coefficients_centroparietal <- data.table::fread("aggregated_data/multitemp_pre200_post800_centroparietal_cluster_omission_trials.csv", quote="") # read data again
data_grand_averages <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_200800_anonym.csv", quote="") # read data again

segment_start <- 200
segment_end <- 800

segment_length <- 250

# axis labels
axis_start <- 200
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- 0
axis_distance <- 50
zero_line_sample_points <- 50

ylim <- c(4,-2) # limits of y-axis in all plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/3)
x_lab_text <- 200

# Open a plot device:
png(filename = paste0("plots/multitemp_results_omission_trials_all_fixed_effects_",as.character(Sys.Date()),".png"),
    width=8, height=7, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(3,2), mai=c(.43,.5,.4,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

#### Plot 1 (topleft): Effects at Frontocentral Cluster | PE ####

# get/rename data
pe_coefficients <- pe_coefficients_frontocentral
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot

## Plot significance indicators

# to loop through fixed effects:
effects <- c("significanteffect_pe", "significanteffect_valence", "significanteffect_interaction")

# Define the mixed color
mid_color <- scales::alpha("#006780", 0.7)
col <- c("purple", mid_color, "cornflowerblue")
pch <- c(15, 15, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}


## Add indicators for model convergence and singular fit

# add yellow circles for singular fit models
if (length(singular)>0){
  for (i in 1:length(singular)){
    points(singular[i], (ylim[1]-(.0375*abs(diff(ylim)))), pch=19, col="coral", xpd=T)
  }
}

# add red circles for non-converged model (overwrite singular fit if it is the same)
if (length(non_convergence)>0){
  for (i in 1:length(non_convergence)){
    points(non_convergence[i], ylim[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
  }
}
## Draw lines for coefficients

lines(pe_coefficients$time, pe_coefficients$coef_valence*2, col=mid_color, lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_valence*2 + pe_coefficients$se_valence, 
                       rev(pe_coefficients$coef_valence*2 - pe_coefficients$se_valence))), col=yarrr::transparent(mid_color, trans.val = .7),xpd=T, border=NA)


lines(pe_coefficients$time, pe_coefficients$coef_interaction_pe_valence, col="cornflowerblue", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_interaction_pe_valence + pe_coefficients$se_interaction_pe_valence, 
                       rev(pe_coefficients$coef_interaction_pe_valence - pe_coefficients$se_interaction_pe_valence))),
        col=yarrr::transparent("cornflowerblue", trans.val = .7),xpd=T, border=NA)


# 
# lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
# lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_pe, col="purple", lwd = 3, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
#                        rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
#                        rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)
# 
# 
polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                       rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent("purple", trans.val = .7),xpd=T, border=NA)


## Annotations
# Axes

axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# y-axis label
title(ylab = expression(beta~"-coefficient"))

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset



#### Plot 2 (middleleft): Grand Averages Separately for Valence ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data
data <- subset(data, omission==1) # only trials with omitted feedback

# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)
# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$valence # separate lines according to which variable?

# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- subset(data_avg, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp


# standard error
temp <- subset(data_se, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


short_pos <- cbind.data.frame(short_pos, short_neg)

short_pos_se <- cbind.data.frame(short_pos_se, short_neg_se)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols=c("blue", "green3"), 
      lwd=3, 
      colsterror=c(yarrr::transparent("blue", trans.val = .7), yarrr::transparent("green3", trans.val = .7))) # topleft
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 3 (bottomleft): Grand Averages Separately for PEs ####

## get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data
data <- subset(data, omission==1) # subset to trials with omitted feedback

## define custom function

x <- seq(-segment_start, segment_end, length.out=segment_length) # define x-axis
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 
data$uniquecond <- data$pe_c_cat # rename to uniquecond to recyclue previous script

## choose electrodes and pipe general info
electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

## aggregate data 

data_avg_all <- data.frame()

# average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- data_avg
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- data_se
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)



gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols= viridis::viridis_pal()(ncol(short_pos)), 
      lwd=3, 
      colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 4 (topright): Effects at Centroparietal Cluster | PE ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot

## Plot significance indicators

# to loop through fixed effects:
effects <- c("significanteffect_pe", "significanteffect_valence", "significanteffect_interaction")

# Define the mixed color
mid_color <- scales::alpha("#006780", 0.7)
col <- c("purple", mid_color, "cornflowerblue")
pch <- c(15, 15, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}

## Add indicators for model convergence and singular fit

# add yellow circles for singular fit models
if (length(singular)>0){
  for (i in 1:length(singular)){
    points(singular[i], (ylim[1]-(.0375*abs(diff(ylim)))), pch=19, col="coral", xpd=T)
  }
}

# add red circles for non-converged model (overwrite singular fit if it is the same)
if (length(non_convergence)>0){
  for (i in 1:length(non_convergence)){
    points(non_convergence[i], ylim[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
  }
}


## Draw lines for coefficients

lines(pe_coefficients$time, pe_coefficients$coef_valence*2, col=mid_color, lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_valence*2 + pe_coefficients$se_valence, 
                       rev(pe_coefficients$coef_valence*2 - pe_coefficients$se_valence))), col=yarrr::transparent(mid_color, trans.val = .7),xpd=T, border=NA)


lines(pe_coefficients$time, pe_coefficients$coef_interaction_pe_valence, col="cornflowerblue", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_interaction_pe_valence + pe_coefficients$se_interaction_pe_valence, 
                       rev(pe_coefficients$coef_interaction_pe_valence - pe_coefficients$se_interaction_pe_valence))),
        col=yarrr::transparent("cornflowerblue", trans.val = .7),xpd=T, border=NA)


# 
# lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
# lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_pe, col="purple", lwd = 3, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
#                        rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
#                        rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)
# 
# 
polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                       rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent("purple", trans.val = .7),xpd=T, border=NA)


## Annotations
# Axes
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset


#### Plot 5 (middleleft): Grand Averages Separately for Valence ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table

aggr_data <- data
#### Plot 2 & 3 (left middle and bottom): Grand Averages | Omitted | PE | Open plot & set global parameters 


# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)


# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


#### Prepare Grand Averages for Frontocentral Electrodes

#data <- aggr_data
data <- subset(data, omission==1) # only trials with omitted feedback

##### Categorize continuous pe variable for separate lines

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$valence # separate lines according to which variable?
#within_plot_variable <- data$valence # separate plots according to which variable?



# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- subset(data_avg, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp


# standard error
temp <- subset(data_se, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


short_pos <- cbind.data.frame(short_pos, short_neg)

short_pos_se <- cbind.data.frame(short_pos_se, short_neg_se)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols=c("blue", "green3"), 
      lwd=3, 
      colsterror=c(yarrr::transparent("blue", trans.val = .7), yarrr::transparent("green3", trans.val = .7))) # topleft
#title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 6 (bottomleft): Grand Averages Separately for PEs ####

## get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data
data <- subset(data, omission==1) # subset to trials with omitted feedback

## define custom function

x <- seq(-segment_start, segment_end, length.out=segment_length) # define x-axis
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 
data$uniquecond <- data$pe_c_cat # rename to uniquecond to recyclue previous script

## choose electrodes and pipe general info
electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

## aggregate data 

data_avg_all <- data.frame()

# average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- data_avg
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- data_se
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)



gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols= viridis::viridis_pal()(ncol(short_pos)), 
      lwd=3, 
      colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
#title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Close plot ####
dev.off() # close device





#### Display | All Fixed Effect Coefficients| -200 to 800 ms || Load data, set general settings and open plot ####

pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre200_post800_frontocentral_cluster_display_trials.csv", quote="") # read data again
pe_coefficients_centroparietal <- data.table::fread("aggregated_data/multitemp_pre200_post800_centroparietal_cluster_display_trials.csv", quote="") # read data again
data_grand_averages <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_only_immediate_downsampled.csv", quote="") # read data again

segment_start <- 200
segment_end <- 800

segment_length <- 250

# axis labels
axis_start <- 200
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- 0
axis_distance <- 50
zero_line_sample_points <- 50

ylim <- c(16,-4) # limits of y-axis in all plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/5)
x_lab_text <- 150

# Open a plot device:
png(filename = paste0("plots/multitemp_results_display_trials_all_fixed_effects_",as.character(Sys.Date()),".png"),
    width=8, height=7, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(3,2), mai=c(.43,.5,.4,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )


#### Plot 1 (topleft): Fixed Effect Coefficients (FC Cluster) ####

# get/rename data
pe_coefficients <- pe_coefficients_frontocentral
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot

## Plot significance indicators

# to loop through fixed effects:
effects <- c("significanteffect_pe", "significanteffect_valence", "significanteffect_interaction")

# Define the mixed color
mid_color <- scales::alpha("#006780", 0.7)
col <- c("purple", mid_color, "cornflowerblue")
pch <- c(15, 15, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}


## Add indicators for model convergence and singular fit

# add yellow circles for singular fit models
if (length(singular)>0){
  for (i in 1:length(singular)){
    points(singular[i], (ylim[1]-(.0375*abs(diff(ylim)))), pch=19, col="coral", xpd=T)
  }
}

# add red circles for non-converged model (overwrite singular fit if it is the same)
if (length(non_convergence)>0){
  for (i in 1:length(non_convergence)){
    points(non_convergence[i], ylim[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
  }
}

## Draw lines for coefficients

lines(pe_coefficients$time, pe_coefficients$coef_valence*2, col=mid_color, lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_valence*2 + pe_coefficients$se_valence, 
                       rev(pe_coefficients$coef_valence*2 - pe_coefficients$se_valence))), col=yarrr::transparent(mid_color, trans.val = .7),xpd=T, border=NA)


lines(pe_coefficients$time, pe_coefficients$coef_interaction_pe_valence, col="cornflowerblue", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_interaction_pe_valence + pe_coefficients$se_interaction_pe_valence, 
                       rev(pe_coefficients$coef_interaction_pe_valence - pe_coefficients$se_interaction_pe_valence))),
        col=yarrr::transparent("cornflowerblue", trans.val = .7),xpd=T, border=NA)


# 
# lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
# lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_pe, col="purple", lwd = 3, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
#                        rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
#                        rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)
# 
# 
polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                       rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent("purple", trans.val = .7),xpd=T, border=NA)


## Annotations

# Axes

axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
#text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# y-axis label
title(ylab = expression(beta~"-coefficient"))

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset



#### Plot 2 (middleleft): Grand Averages Separately for Valence ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table

aggr_data <- data
#### Plot 2 & 3 (left middle and bottom): Grand Averages | Omitted | PE | Open plot & set global parameters 


# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)


# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


#### Prepare Grand Averages for Frontocentral Electrodes

#data <- aggr_data
data <- subset(data, omission==-1) # only trials with presented feedback

##### Categorize continuous pe variable for separate lines

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$valence # separate lines according to which variable?
#within_plot_variable <- data$valence # separate plots according to which variable?



# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- subset(data_avg, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp


# standard error
temp <- subset(data_se, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


short_pos <- cbind.data.frame(short_pos, short_neg)

short_pos_se <- cbind.data.frame(short_pos_se, short_neg_se)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols=c("blue", "green3"), 
      lwd=3, 
      colsterror=c(yarrr::transparent("blue", trans.val = .7), yarrr::transparent("green3", trans.val = .7))) # topleft
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 3 (bottomleft): Grand Averages Separately for PEs ####

## get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data
data <- subset(data, omission==-1) # subset to trials with presented feedback

## define custom function

x <- seq(-segment_start, segment_end, length.out=segment_length) # define x-axis
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 
data$uniquecond <- data$pe_c_cat # rename to uniquecond to recyclue previous script

## choose electrodes and pipe general info
electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

## aggregate data 

data_avg_all <- data.frame()

# average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- data_avg
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- data_se
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)



gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols= viridis::viridis_pal()(ncol(short_pos)), 
      lwd=3, 
      colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 4 (topright): Effects at Centroparietal Cluster | PE ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
plot(1, type="l", ylim=ylim, xlim=c(0,segment_length), axes=F, ylab = NA, xlab = "")#"time, ms") # empty plot

## Plot significance indicators

# to loop through fixed effects:
effects <- c("significanteffect_pe", "significanteffect_valence", "significanteffect_interaction")

# Define the mixed color
mid_color <- scales::alpha("#006780", 0.7)
col <- c("purple", mid_color, "cornflowerblue")
pch <- c(15, 15, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.21*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}

## Add indicators for model convergence and singular fit

# add yellow circles for singular fit models
if (length(singular)>0){
  for (i in 1:length(singular)){
    points(singular[i], (ylim[1]-(.0375*abs(diff(ylim)))), pch=19, col="coral", xpd=T)
  }
}

# add red circles for non-converged model (overwrite singular fit if it is the same)
if (length(non_convergence)>0){
  for (i in 1:length(non_convergence)){
    points(non_convergence[i], ylim[1], pch=19, col="darkgoldenrod1", xpd=T)#col=10)
  }
}


## Draw lines for coefficients

lines(pe_coefficients$time, pe_coefficients$coef_valence*2, col=mid_color, lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_valence*2 + pe_coefficients$se_valence, 
                       rev(pe_coefficients$coef_valence*2 - pe_coefficients$se_valence))), col=yarrr::transparent(mid_color, trans.val = .7),xpd=T, border=NA)


lines(pe_coefficients$time, pe_coefficients$coef_interaction_pe_valence, col="cornflowerblue", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_interaction_pe_valence + pe_coefficients$se_interaction_pe_valence, 
                       rev(pe_coefficients$coef_interaction_pe_valence - pe_coefficients$se_interaction_pe_valence))),
        col=yarrr::transparent("cornflowerblue", trans.val = .7),xpd=T, border=NA)


# 
# lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
# lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_pe, col="purple", lwd = 3, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
#                        rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)
# 
# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
#                        rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)
# 
# 
polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pe + pe_coefficients$se_pe, 
                       rev(pe_coefficients$coef_pe - pe_coefficients$se_pe))), col=yarrr::transparent("purple", trans.val = .7),xpd=T, border=NA)


## Annotations

# Axes
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset


#### Plot 5 (middleleft): Grand Averages Separately for Valence ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data

#### Plot 2 & 3 (left middle and bottom): Grand Averages | Omitted | PE | Open plot & set global parameters 


# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)


# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


#### Prepare Grand Averages for Frontocentral Electrodes

data <- subset(data, omission==-1) # only trials with presented feedback

##### Categorize continuous pe variable for separate lines

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$valence # separate lines according to which variable?

# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- line_variable
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- subset(data_avg, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, uniquecond==1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp


# standard error
temp <- subset(data_se, uniquecond==-1)
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


short_pos <- cbind.data.frame(short_pos, short_neg)

short_pos_se <- cbind.data.frame(short_pos_se, short_neg_se)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols=c("blue", "green3"), 
      lwd=3, 
      colsterror=c(yarrr::transparent("blue", trans.val = .7), yarrr::transparent("green3", trans.val = .7))) # topleft
#title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 6 (bottomleft): Grand Averages Separately for PEs ####

## get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table
aggr_data <- data
data <- subset(data, omission==-1) # subset to trials with presented feedback

## define custom function

x <- seq(-segment_start, segment_end, length.out=segment_length) # define x-axis
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}

## categorize continuous pe variable for separate lines
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 
data$uniquecond <- data$pe_c_cat # rename to uniquecond to recyclue previous script

## choose electrodes and pipe general info
electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

## aggregate data 

data_avg_all <- data.frame()

# average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# positive

# average
temp <- data_avg
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- data_se
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)



gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]),
      sterror = c(short_pos_se[,1:ncol(short_pos)]),
      cols= viridis::viridis_pal()(ncol(short_pos)), 
      lwd=3, 
      colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
#title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Close plot ####
dev.off() # close device





#### Omission | Simple Slopes of PE | -700 to 1500 ms || Load data, set general settings and open plot ####

pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre700_post1500_frontocentral_cluster_omission_trials.csv", quote="") # read data again
pe_coefficients_centroparietal <- data.table::fread("aggregated_data/multitemp_pre700_post1500_centroparietal_cluster_omission_trials.csv", quote="") # read data again
data_grand_averages <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_7001500_anonym.csv", quote="") # read data again

segment_start <- 700
segment_end <- 1500

segment_length <- 550

# axis labels
axis_start <- 500
axis_label_distance <- 500

# axis ticks (in samplepoints)
axis_start_sample_points <- 50
axis_distance <- 125
zero_line_sample_points <- 175

ylim <- c(4,-4) # limits of y-axis in all plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/4)
x_lab_text <- 275


# Open a plot device:
png(filename = paste0("plots/multitemp_results_omission_trials_pe_extended_segment_",as.character(Sys.Date()),".png"),
    width=8, height=7, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(3,2), mai=c(.43,.5,.4,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

#### Plot 1 (topleft): Effects at Frontocentral Cluster | PE ####

# get/rename data
pe_coefficients <- pe_coefficients_frontocentral
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
mid_color <- scales::alpha("#006780", 0.7)
col <- c("blue", "green3", "cyan1")

pch <- c(15, 15, 15)
#pch <- c(23, 23,23,  19, 19, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.28*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}




## Draw lines for coefficients
lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
#lines(pe_coefficients$time, pe_coefficients$coeff_diff, col="cyan1", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
                       rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
                       rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)


# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coeff_diff + pe_coefficients$se_diff, 
#                        rev(pe_coefficients$coeff_diff - pe_coefficients$se_diff))), col=yarrr::transparent("cyan1", trans.val = .7),xpd=T, border=NA)

## Annotations

# Axes
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
#axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance)) # overwrite x-axis (caution!)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# y-axis label
title(ylab = expression(beta~"-coefficient of PE"))

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset



#### Read and prepare data for Grand Averages ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table

aggr_data <- data
#### Plot 2 & 3 (left middle and bottom): Grand Averages | Omitted | PE | Open plot & set global parameters ####


# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)


# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8

  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    #lines(x, aoc, col=cols[2], ...)
    #lines(x, avc, col=cols[3], ...)
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
  
  ## custom y-axis label (mikroV)
 # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


#### Prepare Grand Averages for Frontocentral Electrodes ####
## Prepare grand averages

data <- aggr_data
data <- subset(data, omission==1) # only trials with omitted feedback

##### Categorize continuous pe variable for separate lines

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$pe_c_cat # separate lines according to which variable?
within_plot_variable <- data$valence # separate plots within one png according to which variable?



# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- paste(line_variable, within_plot_variable)
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# Prepare data for plotting

# coding in data$ valence: -1 == negative outcome; 1 == positive outcome

# subset

# positive

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]), sterror = c(short_pos_se[,1:ncol(short_pos)]),  cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
title(ylab = expression(bold("Positive Outcome")), line=2.6) # add y-axis label //expression(bold("Marginal Effect of PE (a.u.)"))
title(ylab = expression("amplitude, "~mu*V), line=1.2)
gravg(x, linecat= c(short_neg[,1:ncol(short_pos)]), sterror= c(short_neg_se[,1:ncol(short_pos)]), cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # bottomleft
title(ylab = expression(bold("Negative Outcome")), line=2.6)#1.5) # add y-axis label
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 4 (topright): Effects at Centroparietal Cluster | PE ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
mid_color <- scales::alpha("#006780", 0.7)
col <- c("blue", "green3", "cyan1")

pch <- c(15, 15, 15)
#pch <- c(23, 23,23,  19, 19, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.28*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}



## Draw lines for coefficients
lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
#lines(pe_coefficients$time, pe_coefficients$coeff_diff, col="cyan1", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
                       rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
                       rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)


# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coeff_diff + pe_coefficients$se_diff, 
#                        rev(pe_coefficients$coeff_diff - pe_coefficients$se_diff))), col=yarrr::transparent("cyan1", trans.val = .7),xpd=T, border=NA)

## Annotations

# Axes
#axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance)) # overwrite x-axis (caution!)
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)

axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset


#### Plot 5 & 6 (right middle and bottom): Prepare Grand Averages for Centroparietal Electrodes ####
## Prepare grand averages

# read data if it is not in the environment
data <- aggr_data
data <- subset(data, omission==1) # only trials with presented feedback

##### Categorize continuous pe variable for separate lines 

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")
#electrodes <- c("F3", "Fz", "F4", "FC1", "FC2", "C3", "Cz", "C4", "P3", "Pz", "P4", "FCz", "CP1", "CP2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?


line_variable <- data$pe_c_cat # separate lines according to which variable?
within_plot_variable <- data$valence # separate plots within one png according to which variable?


# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- paste(line_variable, within_plot_variable)
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# Prepare data for plotting

# coding in data$ valence: -1 == negative outcome; 1 == positive outcome

# subset

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)

gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]), sterror = c(short_pos_se[,1:ncol(short_pos)]),  cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
gravg(x, linecat= c(short_neg[,1:ncol(short_pos)]), sterror= c(short_neg_se[,1:ncol(short_pos)]), cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # bottomleft

#### Close plot ####
dev.off() # close device





#### Clear workspace ####
rm(list=ls())
#### Omission | Simple Slopes of PE | -200 to 800 ms || Load data, set general settings and open plot ####

pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre200_post800_frontocentral_cluster_omission_trials.csv", quote="") # read data again
pe_coefficients_centroparietal <- data.table::fread("aggregated_data/multitemp_pre200_post800_centroparietal_cluster_omission_trials.csv", quote="") # read data again
data_grand_averages <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_200800_anonym.csv", quote="") # read data again

segment_start <- 200
segment_end <- 800

segment_length <- 250

# axis labels
axis_start <- 200
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- 0
axis_distance <- 50
zero_line_sample_points <- 50

ylim <- c(4,-2) # limits of y-axis in all plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/3)
x_lab_text <- 125

# Open a plot device:
png(filename = paste0("plots/multitemp_results_omission_trials_pe_",as.character(Sys.Date()),".png"),
    width=8, height=7, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(3,2), mai=c(.43,.5,.4,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )


#### Plot 1 (topleft): Effects at Frontocentral Cluster | PE ####

# get/rename data
pe_coefficients <- pe_coefficients_frontocentral
pe_coefficients$time <- 1:nrow(pe_coefficients) 
## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
mid_color <- scales::alpha("#006780", 0.7)
col <- c("blue", "green3", "cyan1")

pch <- c(15, 15, 15)
#pch <- c(23, 23,23,  19, 19, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.28*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}



## Draw lines for coefficients
lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
#lines(pe_coefficients$time, pe_coefficients$coeff_diff, col="cyan1", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
                       rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
                       rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)


# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coeff_diff + pe_coefficients$se_diff, 
#                        rev(pe_coefficients$coeff_diff - pe_coefficients$se_diff))), col=yarrr::transparent("cyan1", trans.val = .7),xpd=T, border=NA)
 

## Annotations

# Axes

#axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance)) # overwrite x-axis (caution!)
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)


axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# y-axis label
title(ylab = expression(beta~"-coefficient of PE"))

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset



#### Read and prepare data for Grand Averages ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table

aggr_data <- data
#### Plot 2 & 3 (left middle and bottom): Grand Averages | Omitted | PE | Open plot & set global parameters ####


# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)


# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    #lines(x, aoc, col=cols[2], ...)
    #lines(x, avc, col=cols[3], ...)
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


#### Prepare Grand Averages for Frontocentral Electrodes ####

#data <- aggr_data
data <- subset(data, omission==1) # only trials with omitted feedback

##### Categorize continuous pe variable for separate lines

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$pe_c_cat # separate lines according to which variable?
within_plot_variable <- data$valence # separate plots according to which variable?



# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- paste(line_variable, within_plot_variable)
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# Prepare data for plotting

# coding in data$ valence: -1 == negative outcome; 1 == positive outcome

# subset

# positive

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]), sterror = c(short_pos_se[,1:ncol(short_pos)]),  cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
title(ylab = expression(bold("Positive Outcome")), line=2.6) # add y-axis label //expression(bold("Marginal Effect of PE (a.u.)"))
title(ylab = expression("amplitude, "~mu*V), line=1.2)
gravg(x, linecat= c(short_neg[,1:ncol(short_pos)]), sterror= c(short_neg_se[,1:ncol(short_pos)]), cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # bottomleft
title(ylab = expression(bold("Negative Outcome")), line=2.6)#1.5) # add y-axis label
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 4 (topright): Effects at Centroparietal Cluster | PE ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
mid_color <- scales::alpha("#006780", 0.7)
col <- c("blue", "green3", "cyan1")

pch <- c(15, 15, 15)
#pch <- c(23, 23,23,  19, 19, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.28*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}





## Draw lines for coefficients
lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
#lines(pe_coefficients$time, pe_coefficients$coeff_diff, col="cyan1", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
                       rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
                       rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)


# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coeff_diff + pe_coefficients$se_diff, 
#                        rev(pe_coefficients$coeff_diff - pe_coefficients$se_diff))), col=yarrr::transparent("cyan1", trans.val = .7),xpd=T, border=NA)

## Annotations

# Axes
#axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance)) # overwrite x-axis (caution!)
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)

axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset


#### Plot 5 & 6 (right middle and bottom): Prepare Grand Averages for Centroparietal Electrodes ####
## Prepare grand averages

# read data if it is not in the environment
data <- aggr_data
data <- data.table::setDT(data)
data <- subset(data, omission==1) # only trials with omitted feedback

##### Categorize continuous pe variable for separate lines 

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")
#electrodes <- c("F3", "Fz", "F4", "FC1", "FC2", "C3", "Cz", "C4", "P3", "Pz", "P4", "FCz", "CP1", "CP2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?


line_variable <- data$pe_c_cat # separate lines according to which variable?
within_plot_variable <- data$valence # separate plots according to which variable?


# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- paste(line_variable, within_plot_variable)
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# Prepare data for plotting

# coding in data$ valence: -1 == negative outcome; 1 == positive outcome

# subset

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)

gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]), sterror = c(short_pos_se[,1:ncol(short_pos)]),  cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
gravg(x, linecat= c(short_neg[,1:ncol(short_pos)]), sterror= c(short_neg_se[,1:ncol(short_pos)]), cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # bottomleft


#### Close plot ####
dev.off() # close device





#### Display | Simple Slopes of PE | -200 to 800 ms || Load data, set general settings and open plot ####

pe_coefficients_frontocentral <- data.table::fread("aggregated_data/multitemp_pre200_post800_frontocentral_cluster_display_trials.csv", quote="") # read data again
pe_coefficients_centroparietal <- data.table::fread("aggregated_data/multitemp_pre200_post800_centroparietal_cluster_display_trials.csv", quote="") # read data again
data_grand_averages <- data.table::fread("aggregated_data/eeg_and_pe_data_merged_200800_anonym.csv", quote="") # read data again

segment_start <- 200
segment_end <- 800

segment_length <- 250

# axis labels
axis_start <- 200
axis_label_distance <- 200

# axis ticks (in samplepoints)
axis_start_sample_points <- 0
axis_distance <- 50
zero_line_sample_points <- 50

ylim <- c(16,-4) # limits of y-axis in all plots
ylim_max <- ylim[1]
ylim_min <-  ylim[2]
y_axis_ticks_distance <- abs(diff(ylim)/5)
x_lab_text <- 150

# Open a plot device:
png(filename = paste0("plots/multitemp_results_display_trials_pe_",as.character(Sys.Date()),".png"),
    width=8, height=7, unit="in", res=400)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(3,2), mai=c(.43,.5,.4,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )


#### Plot 1 (topleft): Effects at Frontocentral Cluster | PE ####

# get/rename data
pe_coefficients <- pe_coefficients_frontocentral
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
mid_color <- scales::alpha("#006780", 0.7)
col <- c("blue", "green3", "cyan1")

pch <- c(15, 15, 15)
#pch <- c(23, 23,23,  19, 19, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.28*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}




## Draw lines for coefficients
lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
#lines(pe_coefficients$time, pe_coefficients$coeff_diff, col="cyan1", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
                       rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
                       rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)


# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coeff_diff + pe_coefficients$se_diff, 
#                        rev(pe_coefficients$coeff_diff - pe_coefficients$se_diff))), col=yarrr::transparent("cyan1", trans.val = .7),xpd=T, border=NA)

## Annotations

# Axes

axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)
axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)
# y-axis label
title(ylab = expression(beta~"-coefficient of PE"))

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset

#### Read and prepare data for Grand Averages ####

# get data
data <- data_grand_averages
data <- data.table::setDT(data) # convert to data.table

aggr_data <- data
#### Plot 2 & 3 (left middle and bottom): Grand Averages | Omitted | PE | Open plot & set global parameters ####


# Define x-axis
x <- seq(-segment_start, segment_end, length.out=segment_length)


# Define grand average function for the plot
gravg <- function(x, linecat, sterror, ylim=c(ylim_max,ylim_min), cols=c(1,2), colsterror=c(1,2),
                  xlab="", ylab=expression(mu*V), ...) {
  
  ## set up an empty plot canvas
  plot(1, type="l", col=NA, xlim=range(x), ylim=ylim, axes=FALSE,
       xlab=xlab, ylab=NA, cex.axis=1) # vorher cex.axis=.8
  
  axis(1, line=.7) # x-axis
  abline(h=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  abline(v=0, lty=3) # h=0: horizontal line at 0, lty=3: dashed
  
  
  ## DRAW DATA
  for (i in 1:length(linecat)) {
    # polygon(x)
    polygon(x=c(x, rev(x)), y=c(as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), 
                                rev(as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]))), col=colsterror[i], border=NA, xpd=T)
    #rect(x-0.5, as.numeric(linecat[[i]]) - as.numeric(sterror[[i]]), x+0.5, as.numeric(linecat[[i]]) + as.numeric(sterror[[i]]), col=colsterror[i], border=NA)
  }
  
  for (i in 1:length(linecat)) {
    lines(x, linecat[[i]], col=cols[i])
    
  }
  
  
  
  ## axis with labels, pos=0: at origin, las=2: horizontal labels
  axis(2, at=seq(ylim[2],ylim[1], y_axis_ticks_distance), pos=NA, las=2, cex.axis=1)# cex.axis=.8) # vorher pos=0
  ## small intermediate axis ticks
  axis(2, at=seq(ylim[2],ylim[1], (y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)#cex.axis=.8)
  
  text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)
  
  
  ## custom y-axis label (mikroV)
  # text(x=-650, y=(ylim[2]), expression(mu*V)) # vorher -150
}


#### Prepare Grand Averages for Frontocentral Electrodes ####

#data <- aggr_data
data <- subset(data, omission==-1) # only trials with presented feedback

##### Categorize continuous pe variable for separate lines

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c("F3", "Fz", "F4", "FC1", "FC2")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?

line_variable <- data$pe_c_cat # separate lines according to which variable?
within_plot_variable <- data$valence # separate plots according to which variable?



# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- paste(line_variable, within_plot_variable)
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# Prepare data for plotting

# coding in data$ valence: -1 == negative outcome; 1 == positive outcome

# subset

# positive

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)


gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]), sterror = c(short_pos_se[,1:ncol(short_pos)]),  cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
title(ylab = expression(bold("Positive Outcome")), line=2.6) # add y-axis label //expression(bold("Marginal Effect of PE (a.u.)"))
title(ylab = expression("amplitude, "~mu*V), line=1.2)
gravg(x, linecat= c(short_neg[,1:ncol(short_pos)]), sterror= c(short_neg_se[,1:ncol(short_pos)]), cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # bottomleft
title(ylab = expression(bold("Negative Outcome")), line=2.6)#1.5) # add y-axis label
title(ylab = expression("amplitude, "~mu*V), line=1.2)



#### Plot 4 (topright): Effects at Centroparietal Cluster | PE ####

# read and prepare data
pe_coefficients <- pe_coefficients_centroparietal
pe_coefficients$time <- 1:nrow(pe_coefficients) 

## Apply alpha correction

significanteffect_pe <- which(pe_coefficients$prob_pe < .001) 
significanteffect_valence <- which(pe_coefficients$prob_valence  < .001)
significanteffect_interaction <- which(pe_coefficients$prob_interaction  < .001)

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
mid_color <- scales::alpha("#006780", 0.7)
col <- c("blue", "green3", "cyan1")

pch <- c(15, 15, 15)
#pch <- c(23, 23,23,  19, 19, 15)

for (i in 1:length(effects)){
  
  current_effect <- get(effects[i])
  
  if (length(current_effect>0)){
    for (j in 1: length(current_effect)){
      points(current_effect[j],
             ylim[2]-(0.28*abs(diff(ylim)))+(i*(.05*abs(diff(ylim)))),
             pch=pch[i], bg=col[i], lwd=1.5, col=col[i], cex=1.2, xpd=T)    }
  }
}



## Draw lines for coefficients
lines(pe_coefficients$time, pe_coefficients$coef_pos, col="blue", lwd = 3, xpd=T)
lines(pe_coefficients$time, pe_coefficients$coef_neg, col="green3", lwd = 3, xpd=T)
#lines(pe_coefficients$time, pe_coefficients$coeff_diff, col="cyan1", lwd = 3, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_pos + pe_coefficients$se_pos, 
                       rev(pe_coefficients$coef_pos - pe_coefficients$se_pos))), col=yarrr::transparent("blue", trans.val = .7),border=NA, xpd=T)

polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
        y=as.numeric(c(pe_coefficients$coef_neg + pe_coefficients$se_neg, 
                       rev(pe_coefficients$coef_neg - pe_coefficients$se_neg))), col=yarrr::transparent("green3", trans.val = .7),xpd=T, border=NA)

# polygon(x=as.numeric(c(pe_coefficients$time, rev(pe_coefficients$time))),
#         y=as.numeric(c(pe_coefficients$coeff_diff + pe_coefficients$se_diff, 
#                        rev(pe_coefficients$coeff_diff - pe_coefficients$se_diff))), col=yarrr::transparent("cyan1", trans.val = .7),xpd=T, border=NA)

## Annotations

# Axes
axis(1, at = seq(axis_start_sample_points,segment_length,axis_distance), labels = seq(-axis_start, segment_end, axis_label_distance), line=.7) # overwrite x-axis (caution!)

#text(x_lab_text,(ylim_max+.3*diff(c(ylim_min, ylim_max))), "time, ms", xpd=T, adj=-.5, cex=1.2)

axis(2, at = seq(ylim[2],ylim[1],y_axis_ticks_distance), labels = seq(ylim[2],ylim[1],y_axis_ticks_distance), las=2)

# small intermediate axis ticks
axis(2, at=seq(ylim[2],ylim[1],(y_axis_ticks_distance/2)), pos=NA, labels=NA, tcl=par("tcl")/2, cex.axis=1)

# Time Zero Line
abline(v=zero_line_sample_points, h=0, lty=3) # dashed zerolines at feedback onset


#### Plot 5 & 6 (right middle and bottom): Prepare Grand Averages for Centroparietal Electrodes ####
## Prepare grand averages

# read data if it is not in the environment
data <- aggr_data
data <- data.table::setDT(data)
data <- subset(data, omission==-1) # only trials with presented feedback

##### Categorize continuous pe variable for separate lines 

# Take absolute values of PEs (and we code valence separately)
data$pe_absolute <- abs(data$PE)
data$pe_c_cat <- cut(data$pe_absolute, breaks = 4, include.lowest = TRUE, labels = FALSE) 

#### Insert dataspecific info

electrodes <- c( "CP1", "CP2", "P3", "Pz", "P4")

# data from which electrodes should be used?
start <- -segment_start # first sample point relative to event?
end <- segment_end # last sample point relative to event?


line_variable <- data$pe_c_cat # separate lines according to which variable?
within_plot_variable <- data$valence # separate plots according to which variable?


# name all variables for which separate plots should be created as specified
# before (no change needed when the three variables in paste() were defined):
data$uniquecond <- paste(line_variable, within_plot_variable)
# because former paste-commanded also combined NAs, exclude NAs
data$uniquecond[grep("NA", data$uniquecond)] <- NA 
# create list of unique conditions outside of data without NAs
unique_cond <- unique(data$uniquecond)[!is.na(unique(data$uniquecond))]


#### Aggregate data 

data_avg_all <- data.frame()

# Average across trials of each electrode separately for each participant and each 'unique condition' 
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

## To get averages which are averaged at first across participants for each
## electrode and then averaged across electrodes

# delete id column
data_avg_all_without_id <- data_avg_all[,-which(names(data_avg_all)=="id")]

# average across particpants
data_avg_plot <- aggregate(. ~ uniquecond + electrode, data= data_avg_all_without_id, mean)

# delete electrode column
data_avg_plot <- data_avg_plot[,-which(names(data_avg_plot)=="electrode")]

# average across electrodes
data_avg_plot <- aggregate(. ~ uniquecond, data= data_avg_plot, mean)


## To get standard errors for the average across participants 

# 1. Average across electrode averages for each condition

# Delete electrode column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="electrode")]

# Create averages for each participant separately for the unique condition combinations
data_avg_all <- aggregate(. ~ id + uniquecond, data= data_avg_all, mean)

# Delete id column
data_avg_all <- data_avg_all[,-which(names(data_avg_all)=="id")]

# 2. Average across particpants
data_avg <- c()
data_avg <- aggregate(. ~ uniquecond, data= data_avg_all, mean)

# Define function for standard error computation
std <- function(x) sd(x)/sqrt(length(x))

# Calculate SE for average across participants
data_se <- c()
data_se <- aggregate(. ~ uniquecond, data= data_avg_all, std)


data_avg <- data_avg_plot

# Prepare data for plotting

# coding in data$ valence: -1 == negative outcome; 1 == positive outcome

# subset

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_pos <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))==" 1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_pos_se <- do.call(cbind.data.frame, temp)


# negative

# average
temp <- subset(data_avg, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric) # convert columns to numeric
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0) # apply rolling mean
temp <- do.call(cbind.data.frame, temp) # convert to data frame
short_neg <- temp

# standard error
temp <- subset(data_se, substr(uniquecond, nchar(uniquecond)-1, nchar(uniquecond))=="-1")
temp <- as.data.frame(t(temp)) # transpose (pe categories in columns)
names(temp) <- temp[1,] # first row as column names
temp <- as.data.frame(temp[-1,]) # delete first row and rename data
temp[,1:ncol(temp)] <- sapply(temp[,1:ncol(temp)],as.numeric)
temp <- data.table::frollapply(temp, 10, mean, align="center", fill=0)
short_neg_se <- do.call(cbind.data.frame, temp)

gravg(x, linecat= c(short_pos[,1:ncol(short_pos)]), sterror = c(short_pos_se[,1:ncol(short_pos)]),  cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3, colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # topleft
gravg(x, linecat= c(short_neg[,1:ncol(short_pos)]), sterror= c(short_neg_se[,1:ncol(short_pos)]), cols=viridis::viridis_pal()(ncol(short_pos)), lwd=3,  colsterror=viridis::viridis_pal(alpha=.4)(ncol(short_pos))) # bottomleft

#### Close plot ####
dev.off() # close device





#### Legend ####

png(filename = paste0("plots/legend.png"),
    width=6, height=3, unit="in", res=600)

# Set plot parameters:
# mfcol: multiple plots, mai: plot margins, bottom/left/top/right, tcl: tick mark length
par(mfcol=c(1,1), mai=c(.5,.5,.4,0), mgp=c(1.5,.4,0), tcl=-.25, cex.lab=1.2 )

# second plot: legend
plot(1, type="l", ylim= c(0,10), xlim=c(0,10), axes=F, ylab = NA, xlab = "") # empty plot

# Legend for Simple Slopes of PEs
text(-1, 11.5, as.expression(bquote(bold("Simple Slopes of Absolute PE for"))), pos=4, xpd=T)
legend(-1,11, c("Positive Outcomes", "Negative Outcomes", "Diff. Test of Simple Slopes"),
       col=c("blue","green3", "cyan1"),bty="n", xpd=T, pch=15, #lty=c(1,1),lwd=c(6,6),
       x.intersp=.7)#, cex=.8) # add legend for line colours

# Legend Fixed Effects Coefficients
text(-1, 6.5, as.expression(bquote(bold("Fixed Effects"))), pos=4, xpd=T)
mid_color <- scales::alpha("#006780", 0.7)
legend(-1,6, c(as.expression(bquote("Absolute PE")),
               as.expression(bquote("Feedback Valence")),
               as.expression(bquote("Absolute PE x Feedback Valence")),
               "","Positive Outcomes","Negative Outcomes", "Diff. Test of Simple Slopes"), 
       col=c("purple",mid_color, "cornflowerblue", "white", "blue", "green3", "cyan1"),
       bty="n", xpd=T,
       x.intersp=.7, lty=1, lwd=6)

# Legend for Model diagnostic
text(5, 11.5, expression(bold("Model Diagnostics")), pos=4, xpd=T)
legend(5,11, c("No Convergence", "Singular-fit"),
       col=c("darkgoldenrod1", "coral"), pch=19, bty="n", xpd=T, #title= "Model Diagnostics",
       x.intersp=.7)#, cex=.8) # add legend for line colours

# Legend for PE GAs
text(5, 7, expression(bold("ERPs by Absolute PE")), pos=4, xpd=T)
legend(5, 6.5, c("[      0, <0.25 ]", "[ 0.25, <0.50 ]", "[ 0.50, <0.75 ]", 
                 "[ 0.75,   1      ]"), col=viridis::viridis_pal()(ncol(short_pos)), lty=1, lwd=6, #pch=19, 
       #title=expression(bold("Absolute PE")), 
       bty="n", xpd=T, x.intersp=.7)

# Legend for Valence GAs
text(5, 0.5, as.expression(bquote(bold("ERPs by Feedback Valence"))), pos=4, xpd=T)
legend(5,0, c("Positive Outcome", "Negative Outcome"),
       col=c("blue","green3"),bty="n", xpd=T, lty=c(1,1),lwd=c(6,6),
       x.intersp=.7)#, cex=.8) # add legend for line colours

dev.off() # close device