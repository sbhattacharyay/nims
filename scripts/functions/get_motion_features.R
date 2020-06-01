get_motion_features <- function(temp_var){
  band_powerFeats<-do.call(cbind,temp_var[1:6])
  freq_entropyFeats<-do.call(cbind,temp_var[7:12])
  freq_pairsFeats<-do.call(cbind,temp_var[13:24])
  med_freqFeats<-do.call(cbind,temp_var[25:30])
  smaFeats<-do.call(cbind,temp_var[31:36])
  waveletsFeats<-do.call(cbind,temp_var[37:42])
  
  tempList<-list(band_powerFeats,freq_entropyFeats,freq_pairsFeats,med_freqFeats,smaFeats,waveletsFeats)
  names(tempList)<-c("band_powerFeats","freq_entropyFeats","freq_pairsFeats","med_freqFeats","smaFeats","waveletsFeats")
  
  return(tempList)
}