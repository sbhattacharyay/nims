#### Master Script 2: Multiple Imputation of Missing Accelerometery Values ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Characterize missingness of accelerometry data
# III. Impute totally misssing sensor recordings:
# IV. Perform Amelia II for multiple time-series normal missing data imputation

### I. Initialization
library(tidyverse)
library(forecast)
library(Amelia)
library(imputeTS)
library(parallel)
library(R.matlab)
library(car)
library(bestNormalize)
library(tseries)
library(seastests)
if(.Platform$OS.type == "unix") {
  library(doMC)
} else {
  library(doParallel)
}

# Set the number of parallel cores for parallel tuning
no.parallel.cores <- floor(2 * detectCores() / 3)
if(.Platform$OS.type == "unix") {
  registerDoMC(cores = no.parallel.cores)
} else {
  registerDoParallel(cores = no.parallel.cores)
}

# Load all motion features and compile into one dataframe
if (!exists("all.motion.features")) {
  feature.files <- list.files('../features',pattern=glob2rx("features_*.csv"),full.names = T)
  all.motion.features <- as.data.frame(matrix(nrow = 0, ncol = 12))
  for (curr.feature.file in feature.files) {
    curr.upi.df <- read.csv(curr.feature.file) 
    names(all.motion.features) <- names(curr.upi.df)
    all.motion.features <- rbind(all.motion.features,curr.upi.df)
    print(paste('Feature file',which(feature.files == curr.feature.file),'out of',length(feature.files),'appended'))
  }
}

# Calculate sample size
n <- length(unique(all.motion.features$UPI))

### II. Characterize missingness of accelerometry data:
feature.missingness <- all.motion.features %>% 
  group_by(UPI,Feature) %>%
  summarise(n = n(),
            missingBed = sum(is.na(Bed)),
            missingLA = sum(is.na(LA)),
            missingLE = sum(is.na(LE)),
            missingLW = sum(is.na(LW)),
            missingRA = sum(is.na(RA)),
            missingRE = sum(is.na(RE)),
            missingRW = sum(is.na(RW))) %>%
  rowwise() %>%
  mutate_at(vars(starts_with('missing')),~ .x/n) %>%
  pivot_longer(names_prefix = 'missing',cols = starts_with('missing'),
               values_to = 'perc_missing',names_to='sensor')

# Extract totally missing set
totally.missing.set <- feature.missingness %>% 
  filter(perc_missing == 1) %>%
  select(-c(n,perc_missing))

### III. Impute totally misssing sensor recordings:
# First separate out missing upper extremity indices:
totally.missing.UE <- totally.missing.set %>% filter(sensor%in%c('LE','LW','RE','RW'))
uniq.UE.combos <- unique(totally.missing.UE[,c('sensor','Feature')])
UE.mdls <- vector(mode = "list")
UE.bxcx <- vector(mode = "list")

# Iterate through unique combinations of missing upper extremities and features
for (i in 1:nrow(uniq.UE.combos)){
  if(uniq.UE.combos$sensor[i] %in% c(LE,LW)) {
    filtered.set <- all.motion.features %>% drop_na(LE,LW) %>% filter(Feature == uniq.UE.combos$Feature[i])
    # Train and apply Box Cox transformation to non-missing upper-extremity values
    LE.bxcx <- boxcox(filtered.set$LE,standardize = TRUE)
    LW.bxcx <- boxcox(filtered.set$LW,standardize = TRUE)
    UE.bxcx[[i]]<- list(LE.bxcx,LW.bxcx)
    filtered.set$LE <- LE.bxcx$x.t
    filtered.set$LW <- LW.bxcx$x.t
    if (uniq.UE.combos$sensor[i] == 'LE') {
      UE.mdl <- lm(LE ~ LW, filtered.set)
    } else {
      UE.mdl <- lm(LW ~ LE, filtered.set)
    }
  } else if (uniq.UE.combos$sensor[i] %in% c(RE,RW)) {
    filtered.set <- all.motion.features %>% drop_na(RE,RW) %>% filter(Feature == uniq.UE.combos$Feature[i])
    # Train and apply Box Cox transformation to non-missing upper-extremity values
    RE.bxcx <- boxcox(filtered.set$RE,standardize = TRUE)
    RW.bxcx <- boxcox(filtered.set$RW,standardize = TRUE)
    UE.bxcx[[i]]<- list(RE.bxcx,RW.bxcx)
    filtered.set$RE <- RE.bxcx$x.t
    filtered.set$RW <- RW.bxcx$x.t
    if (uniq.UE.combos$sensor[i] == 'RE') {
      UE.mdl <- lm(RE ~ RW, filtered.set)
    } else {
      UE.mdl <- lm(RW ~ RE, filtered.set)
    }
  }
  UE.mdls[[i]] <- UE.mdl
  print(paste('combination no',i,'complete'))
}

# Impute all totally missing upper extremity values possible with regression models:
for (i in 1:nrow(totally.missing.UE)){
  curr.rows <- which(all.motion.features$UPI == totally.missing.UE$UPI[i] & all.motion.features$Feature == totally.missing.UE$Feature[i])
  curr.df <- all.motion.features %>% filter(UPI == totally.missing.UE$UPI[i] & Feature == totally.missing.UE$Feature[i])
  uniqIdx <- which(uniq.UE.combos$sensor == totally.missing.UE$sensor[i] & uniq.UE.combos$Feature == totally.missing.UE$Feature[i])
  curr.UE.mdl <- UE.mdls[[uniqIdx]]
  curr.E.bxcx <- UE.bxcx[[uniqIdx]][[1]]
  curr.W.bxcx <- UE.bxcx[[uniqIdx]][[2]]
  if (totally.missing.UE$sensor[i] == 'LE'){
    curr.df$LW <- predict(curr.W.bxcx,newdata = curr.df$LW)
    tf.outs <- predict(curr.UE.mdl,newdata = curr.df)
    outs <- predict(curr.E.bxcx,newdata = tf.outs,inverse = TRUE)
  } else if (totally.missing.UE$sensor[i] == 'LW') {
    curr.df$LE <- predict(curr.E.bxcx,newdata = curr.df$LE)
    tf.outs <- predict(curr.UE.mdl,newdata = curr.df)
    outs <- predict(curr.W.bxcx,newdata = tf.outs,inverse = TRUE)    
  } else if (totally.missing.UE$sensor[i] == 'RW') {
    curr.df$RW <- predict(curr.W.bxcx,newdata = curr.df$RW)
    tf.outs <- predict(curr.UE.mdl,newdata = curr.df)
    outs <- predict(curr.E.bxcx,newdata = tf.outs,inverse = TRUE)
  } else if (totally.missing.UE$sensor[i] == 'RW') {
    curr.df$RE <- predict(curr.E.bxcx,newdata = curr.df$RE)
    tf.outs <- predict(curr.UE.mdl,newdata = curr.df)
    outs <- predict(curr.W.bxcx,newdata = tf.outs,inverse = TRUE)   
  }
  all.motion.features[curr.rows,totally.missing.UE$sensor[i]] <- outs
}
rm(UE.mdls,UE.mdl,UE.bxcx,curr.E.bxcx,curr.W.bxcx,curr.UE.mdl,LE.bxcx,LW.bxcx,RE.bxcx,RW.bxcx)

# Impute totally missing bed values by randomly sampling (with replacement) from available bed data of the same sensor:
totally.missing.bed <- totally.missing.set %>% filter(sensor == 'Bed')
for (i in 1:nrow(totally.missing.bed)){
  curr.rows <- which(all.motion.features$UPI == totally.missing.bed$UPI[i] & all.motion.features$Feature == totally.missing.bed$Feature[i])
  filtered.bed.set <- all.motion.features %>% drop_na(Bed) %>% filter(Feature == totally.missing.bed$Feature[i])
  all.motion.features[curr.rows,'Bed'] <- sample(x =filtered.bed.set$Bed,size = length(curr.rows),replace = TRUE)
}
rm(filtered.bed.set)

# Imputation of totally-missing lower extremity recordings
totally.missing.LE <- totally.missing.set %>% filter(sensor%in%c('LA','RA'))
uniq.LE.combos <- unique(totally.missing.LE[,c('sensor','Feature')])
LE.mdls <- vector(mode = "list")
LE.bxcx <- vector(mode = "list")

# Train and apply box-cox filter to normalize dataset
for (i in 1:nrow(uniq.LE.combos)){
  filtered.set <- all.motion.features %>% drop_na(LA,RA) %>% filter(Feature == uniq.LE.combos$Feature[i]) 
  LA.bxcx <- boxcox(filtered.set$LA,standardize = TRUE)
  RA.bxcx <- boxcox(filtered.set$RA,standardize = TRUE)
  LE.bxcx[[i]]<- list(LA.bxcx,RA.bxcx)
  filtered.set$LA <- LA.bxcx$x.t
  filtered.set$RA <- RA.bxcx$x.t
  if (uniq.LE.combos$sensor[i] == 'LA'){
    LEmdl <- lm(LA ~ RA, filtered.set)
  } else {
    LEmdl <- lm(RA ~ LA, filtered.set)
  }
  LE.mdls[[i]] <- LEmdl
  print(paste('combination no',i,'complete'))
}

# Impute totally missing lower extremities
for (i in 1:nrow(totally.missing.LE)){
  curr.rows <- which(all.motion.features$UPI == totally.missing.LE$UPI[i] & all.motion.features$Feature == totally.missing.LE$Feature[i])
  curr.df <- all.motion.features %>% filter(UPI == totally.missing.LE$UPI[i] & Feature == totally.missing.LE$Feature[i])
  uniqIdx <- which(uniq.LE.combos$sensor == totally.missing.LE$sensor[i] & uniq.LE.combos$Feature == totally.missing.LE$Feature[i])
  curr.mdl <- LE.mdls[[uniqIdx]]
  curr.LA.bxcx <- LE.bxcx[[uniqIdx]][[1]]
  curr.RA.bxcx <- LE.bxcx[[uniqIdx]][[2]]
  if (totally.missing.LE$sensor[i] == 'LA'){
    curr.df$RA <- predict(curr.RA.bxcx,newdata = curr.df$RA)
    tf.outs <- predict(curr.mdl,newdata = curr.df)
    outs <- predict(curr.LA.bxcx,newdata = tf.outs,inverse = TRUE) 
  } else if (totally.missing.LE$sensor[i] == 'RA'){
    curr.df$LA <- predict(curr.LA.bxcx,newdata = curr.df$LA)
    tf.outs <- predict(curr.mdl,newdata = curr.df)
    outs <- predict(curr.RA.bxcx,newdata = tf.outs,inverse = TRUE) 
  }
  all.motion.features[curr.rows,totally.missing.LE$sensor[i]] <- outs
}
rm(LE.mdls,LEmdl,LE.bxcx,curr.LA.bxcx,curr.RA.bxcx,curr.mdl,LA.bxcx,RA.bxcx,filtered.set)

### IV. Perform Amelia II for multiple time-series normal missing data imputation
sensor.placements <- c("Bed","RE","LE","RW","LW","RA","LA")
stored.amelias <- vector(mode = "list")
stored.bxcx <- vector(mode = "list")
for (i in 1:length(unique(all.motion.features$Feature))){
  curr.amelia.df <- data.frame(matrix(ncol = ncol(all.motion.features), nrow = 0))
  names(curr.amelia.df) <- names(all.motion.features)
  amelia.bxcx <- vector(mode = "list")
  for (curr.UPI in unique(all.motion.features$UPI)){
    curr.df <- all.motion.features %>% filter(UPI == curr.UPI & Feature == unique(all.motion.features$Feature)[i])
    temp.bxcx <- vector(mode = "list")
    # transform each stream with univariate boxcox
    for (k in sensor.placements) {
      curr.bxcx <- boxcox(curr.df[,k],standardize = TRUE)
      curr.df[,k] <- curr.bxcx$x.t
      temp.bxcx[[k]] <- curr.bxcx
    }
    amelia.bxcx[[curr.UPI]] <- temp.bxcx
    curr.amelia.df <- rbind(curr.amelia.df,curr.df)
  }
  curr.amelia.df <- curr.amelia.df %>% dplyr::select(-Feature)
  curr.amelia <- amelia(curr.amelia.df, m = 9, ts = "RecordingIdx", cs ="UPI",polytime=2,intercs = FALSE, p2s = 2)
  stored.amelias[[i]] <- curr.amelia
  stored.bxcx[[i]] <- amelia.bxcx
  print(paste('Feature no.',i,'complete'))
}

dir.create('../features/01_imputed_features',showWarnings = FALSE)
for (i in 1:length(stored.amelias)){
  curr.amelia <- stored.amelias[[i]]
  curr.bxcx <- stored.bxcx[[i]]
  for(l in 1:curr.amelia$m){
    curr.imp <- curr.amelia$imputations[[l]]
    for (curr.UPI in unique(all.motion.features$UPI)){
      curr.imp.pt <- curr.imp %>% filter(UPI == curr.UPI)
      rows.for.change <- which(curr.imp$UPI == curr.UPI)
      temp.bxcx <- curr.bxcx[[curr.UPI]]
      for (k in sensor.placements){
        curr.vec <- predict(temp.bxcx[[k]],newdata = curr.imp.pt[,k],inverse = TRUE)
        if (sum(is.na(curr.vec)) > 0){
          curr.vec <- na_interpolation(curr.vec, option = "linear")
        }
        curr.imp.pt[,k] <- curr.vec
      }
      curr.imp[rows.for.change,] <- curr.imp.pt
    }
      fileName <- paste0(unique(all.motion.features$Feature)[i],"_",l,".csv")
      write.csv(curr.imp,file.path("../features/01_imputed_features",fileName),row.names = F)
  }
  print(paste('Feature no.',i,'complete'))
}