source("./functions/generateRootDir.R")

viz_lol_2D <- function(Y_LOL,Y,directory,height = 1500, width = 1800, ptSize = 20){

  Y<-Y[!is.na(Y)]  
  rootDir<-generateRootDir(file.path("LOL",directory))
  
  # band_power

  data <- data.frame(x1=Y_LOL$band_powerLOL$Xr[,1], x2=Y_LOL$band_powerLOL$Xr[,2],x3=Y_LOL$band_powerLOL$Xr[,3],y=Y)
  viz<-ggplot(data, aes(x=x1, y=x2, color=y)) +
    geom_point() +
    xlab("Projection 1") +
    ylab("Projection 2") +
    ggtitle(paste("LOL: Band Power",directory)) +
    labs(color = directory)
  fileName <- file.path(rootDir,"band_power.png")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  print(viz)
  dev.off()

  # freq_entropy
  
  data <- data.frame(x1=Y_LOL$freq_entropyLOL$Xr[,1], x2=Y_LOL$freq_entropyLOL$Xr[,2],x3=Y_LOL$freq_entropyLOL$Xr[,3],y=Y)
  viz<-ggplot(data, aes(x=x1, y=x2, color=y)) +
    geom_point() +
    xlab("Projection 1") +
    ylab("Projection 2") +
    ggtitle(paste("LOL: Freq Entropy",directory)) +
    labs(color = directory)
  fileName <- file.path(rootDir,"freq_entropy.png")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  print(viz)
  dev.off()
  
  # freq_pairs
  
  data <- data.frame(x1=Y_LOL$freq_pairsLOL$Xr[,1], x2=Y_LOL$freq_pairsLOL$Xr[,2],x3=Y_LOL$freq_pairsLOL$Xr[,3],y=Y)
  viz<-ggplot(data, aes(x=x1, y=x2, color=y)) +
    geom_point() +
    xlab("Projection 1") +
    ylab("Projection 2") +
    ggtitle(paste("LOL: Freq Pairs",directory)) +
    labs(color = directory)
  fileName <- file.path(rootDir,"freq_pairs.png")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  print(viz)
  dev.off()
  
  # med_freq
  
  data <- data.frame(x1=Y_LOL$med_freqLOL$Xr[,1], x2=Y_LOL$med_freqLOL$Xr[,2],x3=Y_LOL$med_freqLOL$Xr[,3],y=Y)
  viz<-ggplot(data, aes(x=x1, y=x2, color=y)) +
    geom_point() +
    xlab("Projection 1") +
    ylab("Projection 2") +
    ggtitle(paste("LOL: Med Freq",directory)) +
    labs(color = directory)
  fileName <- file.path(rootDir,"med_freq.png")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  print(viz)
  dev.off()
  
  # sma
  
  data <- data.frame(x1=Y_LOL$smaLOL$Xr[,1], x2=Y_LOL$smaLOL$Xr[,2],x3=Y_LOL$smaLOL$Xr[,3],y=Y)
  viz<-ggplot(data, aes(x=x1, y=x2, color=y)) +
    geom_point() +
    xlab("Projection 1") +
    ylab("Projection 2") +
    ggtitle(paste("LOL: SMA",directory)) +
    labs(color = directory)
  fileName <- file.path(rootDir,"sma.png")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  print(viz)
  dev.off()
  
  # wavelets
  
  data <- data.frame(x1=Y_LOL$waveletsLOL$Xr[,1], x2=Y_LOL$waveletsLOL$Xr[,2],x3=Y_LOL$waveletsLOL$Xr[,3],y=Y)
  viz<-ggplot(data, aes(x=x1, y=x2, color=y)) +
    geom_point() +
    xlab("Projection 1") +
    ylab("Projection 2") +
    ggtitle(paste("LOL: Wavelets",directory)) +
    labs(color = directory)
  fileName <- file.path(rootDir,"wavelets.png")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  print(viz)
  dev.off()
}