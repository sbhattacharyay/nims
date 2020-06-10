find_thresholds <- function(df,sma_thresh){
  norm_vec <- function(x) sqrt(sum(x^2))
  
  smaBedProps <- df %>% filter(featureType == "sma") %>% group_by(ptIdx) %>% summarise(propNM = sum(Bed > sma_thresh)/n(),.groups = 'drop')
  
  # We now wish to find a threshold for each of the remaining feature spaces that would minimize the distance to the SMA active bed proportions
  feature_thresholds <- rep(NA,length(unique(df$featureType)))
  
  #Finding Band Power threshold:
  j<-1
  temp_feature_df <- df %>% filter(featureType == "band_power") %>% group_by(ptIdx)
  min_norm <- .Machine$double.xmax
  opt_k = NA;
  for (k in seq(0,0.15,by=.001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(smaBedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Freq Entropy threshold:
  j<-2
  temp_feature_df <- df %>% filter(featureType == "freq_entropy") %>% group_by(ptIdx)
  min_norm <- .Machine$double.xmax
  opt_k = NA;
  for (k in seq(1.5,1.7,by=0.001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed < k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(smaBedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Freq Pairs 1 threshold:
  j<-3
  temp_feature_df <- df %>% filter(featureType == "freq_pairs1") %>% group_by(ptIdx)
  min_norm <- .Machine$double.xmax
  opt_k = NA;
  for (k in seq(0.003,0.02,by=0.001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(smaBedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Freq Pairs 2 threshold:
  j<-4
  temp_feature_df <- df %>% filter(featureType == "freq_pairs2") %>% group_by(ptIdx)
  min_norm <- .Machine$double.xmax
  opt_k = NA;
  for (k in seq(0.001,0.010,by=0.0001)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(smaBedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  #Finding Med Freq threshold:
  j<-5
  temp_feature_df <- df %>% filter(featureType == "med_freq") %>% group_by(ptIdx)
  min_norm <- .Machine$double.xmax
  opt_k = NA;
  for (k in seq(1,3.5,by=0.01)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed < k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(smaBedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;

  feature_thresholds[6] <- sma_thresh;
  
  #Finding Wavelet threshold:
  j<-7
  temp_feature_df <- df %>% filter(featureType == "wavelets") %>% group_by(ptIdx)
  min_norm <- .Machine$double.xmax
  opt_k = NA;
  for (k in seq(.5,10,by=.1)){
    currFeatProps <- temp_feature_df %>% summarise(propNM = sum(Bed > k)/n(),.groups = 'drop')
    curr_norm<-norm_vec(smaBedProps$propNM - currFeatProps$propNM)
    if (curr_norm < min_norm){
      opt_k <- k
      min_norm <- curr_norm;
    }
  }
  feature_thresholds[j] <- opt_k;
  
  return(feature_thresholds)
}