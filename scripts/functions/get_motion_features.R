get_motion_features <- function(temp_var){
  band_powerFeats<-do.call(cbind,temp_var[1:69])
  freq_entropyFeats<-do.call(cbind,temp_var[70:138])
  freq_pairsFeats<-do.call(cbind,temp_var[139:276])
  med_freqFeats<-do.call(cbind,temp_var[277:345])
  smaFeats<-do.call(cbind,temp_var[346:414])
  waveletsFeats<-do.call(cbind,temp_var[415:483])
  
  tempList<-list(band_powerFeats,freq_entropyFeats,freq_pairsFeats,med_freqFeats,smaFeats,waveletsFeats)
  names(tempList)<-c("band_powerFeats","freq_entropyFeats","freq_pairsFeats","med_freqFeats","smaFeats","waveletsFeats")
  
  return(tempList)
}