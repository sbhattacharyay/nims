if (!require(lolR)) install_github('neurodata/lol', build_vignettes=TRUE, force=TRUE)

cv_lol_project_motion_features <- function(motion_features,Y,r,logicalTrain){
  
  logicalIdx <- logicalTrain & !is.na(Y)
  
  band_powerLOL <- lol.project.lol(motion_features$band_powerFeats[logicalIdx,], Y[logicalIdx], r)
  freq_entropyLOL <- lol.project.lol(motion_features$freq_entropyFeats[logicalIdx,], Y[logicalIdx], r)
  freq_pairsLOL <- lol.project.lol(motion_features$freq_pairsFeats[logicalIdx,], Y[logicalIdx], r)
  med_freqLOL <- lol.project.lol(motion_features$med_freqFeats[logicalIdx,], Y[logicalIdx], r)
  smaLOL <- lol.project.lol(motion_features$smaFeats[logicalIdx,], Y[logicalIdx], r)
  waveletsLOL <- lol.project.lol(motion_features$waveletsFeats[logicalIdx,], Y[logicalIdx], r)
  
  tempList<-list(band_powerLOL,freq_entropyLOL,freq_pairsLOL,med_freqLOL,smaLOL,waveletsLOL)
  
  names(tempList)<-c("band_powerLOL","freq_entropyLOL","freq_pairsLOL","med_freqLOL","smaLOL","waveletsLOL")
  
  return(tempList)
}