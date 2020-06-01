if (!require(lolR)) install_github('neurodata/lol', build_vignettes=TRUE, force=TRUE)

lol_project_motion_features <- function(motion_features,Y,r){
  
  band_powerLOL <- lol.project.lol(motion_features$band_powerFeats[!is.na(Y),], Y[!is.na(Y)], r)
  freq_entropyLOL <- lol.project.lol(motion_features$freq_entropyFeats[!is.na(Y),], Y[!is.na(Y)], r)
  freq_pairsLOL <- lol.project.lol(motion_features$freq_pairsFeats[!is.na(Y),], Y[!is.na(Y)], r)
  med_freqLOL <- lol.project.lol(motion_features$med_freqFeats[!is.na(Y),], Y[!is.na(Y)], r)
  smaLOL <- lol.project.lol(motion_features$smaFeats[!is.na(Y),], Y[!is.na(Y)], r)
  waveletsLOL <- lol.project.lol(motion_features$waveletsFeats[!is.na(Y),], Y[!is.na(Y)], r)
  
  tempList<-list(band_powerLOL,freq_entropyLOL,freq_pairsLOL,med_freqLOL,smaLOL,waveletsLOL)
  
  names(tempList)<-c("band_powerLOL","freq_entropyLOL","freq_pairsLOL","med_freqLOL","smaLOL","waveletsLOL")
  
  return(tempList)
}