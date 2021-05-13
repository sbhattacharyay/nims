#### Supplementary function: Find limits of static activity for non-SMA motion 
#### features based on equal distributions to the literature SMA threshold ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk

find_thresholds <- function(df,SMA_thresh){
  norm_vec <- function(x) sqrt(sum(x^2))
  
  SMABedProps <- df %>% filter(Feature == "SMA") %>% group_by(UPI) %>% summarise(propNM = sum(Bed > SMA_thresh)/n(),.groups = 'drop')
  
  # We now wish to find a threshold for each of the remaining feature spaces that would minimize the distance to the SMA active bed proportions
  feature_thresholds <- rep(NA,length(unique(df$Feature)))
  
  #Finding Band Power threshold:
  j<-1
  temp_feature_df <- df %>% filter(Feature == "BPW") %>% group_by(UPI)
  min_norm <- .Machine$double.xmax
  opt_k <- NA
  for (k in seq(0,0.15,by=.001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(SMABedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Freq Entropy threshold:
  j<-2
  temp_feature_df <- df %>% filter(Feature == "FDE") %>% group_by(UPI)
  min_norm <- .Machine$double.xmax
  opt_k <- NA
  for (k in seq(1.5,1.7,by=0.001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed < k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(SMABedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Freq Pairs 1 threshold:
  j<-3
  temp_feature_df <- df %>% filter(Feature == "HLF_l") %>% group_by(UPI)
  min_norm <- .Machine$double.xmax
  opt_k <- NA
  for (k in seq(0.003,0.02,by=0.001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(SMABedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Freq Pairs 2 threshold:
  j<-4
  temp_feature_df <- df %>% filter(Feature == "HLF_h") %>% group_by(UPI)
  min_norm <- .Machine$double.xmax
  opt_k <- NA
  for (k in seq(0.001,0.010,by=0.0001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(SMABedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Med Freq threshold:
  j<-5
  temp_feature_df <- df %>% filter(Feature == "MFR") %>% group_by(UPI)
  min_norm <- .Machine$double.xmax
  opt_k <- NA
  for (k in seq(1,3.5,by=0.01)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed < k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(SMABedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;

  feature_thresholds[6] <- SMA_thresh;
  
  #Finding Wavelet threshold:
  j<-7
  temp_feature_df <- df %>% filter(Feature == "WVL") %>% group_by(UPI)
  min_norm <- .Machine$double.xmax
  opt_k <- NA
  for (k in seq(.5,10,by=.1)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(SMABedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  return(feature_thresholds)
}