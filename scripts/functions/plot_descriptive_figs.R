library(beanplot)
library(readxl)
library(tidyr)
library(dplyr)

source("./functions/generateRootDir.R")
source("./functions/plotBeanplots.R")

plotDescriptiveFigs <- function(patient_clinical_data,clinicalVariableList, exclusions = NULL, height = 3000, width = 3600, ptSize = 20,titled = TRUE, legend = TRUE) {
  
  ### Load the list of variables to be plotted
  
  # Mess-Around Zone:
  continuousVariableDF <- filter(clinicalVariableList, type %in% c("numeric", "integer"))
  
  if (!is.null(exclusions)) {
    continuousVariableDF <- continuousVariableDF[!continuousVariableDF$variableSubtitle %in% exclusions,]
  }
  
  ### Create Graph directory
  rootDir <- generateRootDir("")
  
  ###  histogram of age
  
  fileName <- paste(rootDir,"/age_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  hist(patient_clinical_data$age,
       col = c("darkblue"),
       xlab = "Age",
       main = "Histogram of Age",
       xlim = c(0,100),
       xaxt = 'n')
  axis(1, at = seq(0, 100, by = 5), las=2)

  dev.off()
  
  ###  histogram of icu_duration
  fileName <- paste(rootDir,"/icuDuration_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  
  barplot(table(patient_clinical_data$los),
          col = c("darkblue"),
          ylab = "Frequency",
          xlab = "ICU duration (days)",
          main = "Histogram of ICU duration")
  
  dev.off()
  
  
  ###  histogram of alive or dead at discharge
  fileName <- paste(rootDir,"/death_dis_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  
  plot(patient_clinical_data$death,
       col = c("#A3C1AD", "darkblue"),
       ylab = "Frequency",
       ylim = c(0,60),
       xlab = "Classifier",
       main = "Histogram of outcome: Death at discharge")
  
  dev.off()
  
  ###  histogram of favorable or unfavorable at discharge
  fileName <- paste(rootDir,"/fav_dis_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  
  plot(patient_clinical_data$favorable,
       col = c("#A3C1AD", "darkblue"),
       ylab = "Frequency",
       ylim = c(0,60),
       xlab = "Classifier",
       main = "Histogram of outcome: Favorable at discharge")
  
  dev.off()
  
  ###  histogram of alive or death at 12 months
  fileName <- paste(rootDir,"/death_12mo_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  
  plot(patient_clinical_data$death_12mo,
       col = c("#A3C1AD", "darkblue"),
       ylab = "Frequency",
       ylim = c(0,60),
       xlab = "Classifier",
       main = "Histogram of outcome: Death at 12 months")
  
  dev.off()
  
  ###  histogram of favorable or unfavorable at 12 months
  fileName <- paste(rootDir,"/fav_12mo_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  
  plot(patient_clinical_data$favorable_12mo,
       col = c("#A3C1AD", "darkblue"),
       ylab = "Frequency",
       ylim = c(0,60),
       xlab = "Classifier",
       main = "Histogram of outcome: Favorable at 12 months")
  
  dev.off()
  
  ###  histogram of GOSE at discharge
  fileName <- paste(rootDir,"/gose_dis_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  
  plot(as.factor(patient_clinical_data$gose),
       col = c("#A3C1AD", "darkblue"),
       ylab = "Frequency",
       ylim = c(0,30),
       xlab = "Classifier",
       main = "Histogram of outcome: GOSE at Discharge")
  
  dev.off()
  
  ###  histogram of GOSE at 12 months
  fileName <- paste(rootDir,"/gose_12mo_Hist",
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = ptSize)
  par(mfrow = c(1,1))
  
  plot(as.factor(patient_clinical_data$gose_12mo),
       col = c("#A3C1AD", "darkblue"),
       ylab = "Frequency",
       ylim = c(0,30),
       xlab = "Classifier",
       main = "Histogram of outcome: GOSE at 12 months")
  
  dev.off()
  
  ###  Plot the Beanplots
  
  plotList <- unique(dplyr::select(continuousVariableDF, variableSubtitle, figTitle, type, units))
  rootDir<-generateRootDir("beanplots")
  
  for(i in 1:nrow(plotList)) {
    fileName <- paste(rootDir,"/",
                      plotList$variableSubtitle[i],
                      "_Beanplot",
                      ".png",
                      sep = "")
    png(file = fileName,
        width = width,
        height = height,
        res = 300,
        pointsize = ptSize)
    par(mfrow = c(1,1))
    
    plotBeanplots(variable = plotList$variableSubtitle[i], 
                  plotName = plotList$figTitle[i], 
                  plotUnits = plotList$units[i],
                  type = plotList$type[i], 
                  inputData = patient_clinical_data,
                  legend = legend,
                  titled = titled)
    
    dev.off()
  }
  
  
  ### Plot ventilated data
  aliveID <- dataSet$alive_dead_icu == "A"
  deadID <- dataSet$alive_dead_icu != "A"
  
  #pull out ventilated data into a data frame based on the classifier
  generateVentTable <- function(classifierList) {
    #pull out ventilated data for the classifier
    tempVentData <- data.frame(summary(inputData[classifierList,grepl("Ventilated",colnames(inputData))]))
    
    #reformat dataframe
    tempVentData$Type <- c("No","Yes","NA")
    tempVentData$Freq <- lapply(tempVentData$Freq, gsub, pattern = "Yes|No|:|NA's", replacement = "")
    tempVentData$Freq <- as.integer(lapply(tempVentData$Freq, trimws, which = "both"))
    tempVentData <- tempVentData[,2:length(tempVentData)]
    colnames(tempVentData) <- c("Day","Freq","Type")
    
    #return vetilation table
    return(tempVentData)
  }
  
  #function that returns the percentages of ventilated patients based on the classified data frame
  getClassifierPercs <- function(tempVentData) {
    #return % of ventilated for each day
    return(tempVentData$Freq[tempVentData$Type == "Yes"] / (tempVentData$Freq[tempVentData$Type == "Yes"] + tempVentData$Freq[tempVentData$Type == "No"]))
  }
  
  getVentilatedPercs <- function(tempVentDataAlive, tempVentDataDead, ventilated) {
    tempVentDataAlive$Freq[tempVentDataDead$Type == ventilated] / (tempVentDataDead$Freq[tempVentDataAlive$Type == ventilated] + tempVentDataAlive$Freq[tempVentDataDead$Type == ventilated])
  }
  
  ventData1 = data.frame(list(c(getVentilatedPercs(generateVentTable(aliveID),generateVentTable(!aliveID),"Yes"), getVentilatedPercs(generateVentTable(aliveID),generateVentTable(!aliveID),"No")), 
                              c(replicate(5,"Ventilated"),replicate(5,"Not Ventilated")),
                              c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5")))
  
  colnames(ventData1) = c("Percent_Alive", "Ventilated", "Day")
  
  ventData2 = data.frame(list(c(getClassifierPercs(generateVentTable(aliveID)), getClassifierPercs(generateVentTable(!aliveID))), 
                              c(replicate(5,"Alive"),replicate(5,"Dead")),
                              c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5")))
  
  colnames(ventData2) = c("Percent_Ventilated", "Classifier", "Day")
  
  #draw the graphs
  library(lattice)
  
  fileName <- paste(rootDir,
                    "/Ventilated",
                    "_Barplot_V1_",
                    Sys.Date(),
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = 12)
  par(mfrow = c(1,1))
  
  barchart(Percent_Alive ~ Day, 
           groups=Ventilated, 
           ventData1, 
           auto.key = list(columns = 3),
           par.settings = list(
             superpose.polygon=list(col=c("darkgrey","darkred"), border="black"),
             strip.background=list(col="red"),
             strip.border=list(col="black")),
           ylim = c(0,1),
           main="Percentage dead based on ventilation status")
  
  dev.off()
  
  fileName <- paste(rootDir,
                    "/Ventilated",
                    "_Barplot_V2_",
                    Sys.Date(),
                    ".png",
                    sep = "")
  png(file = fileName,
      width = width,
      height = height,
      res = 300,
      pointsize = 12)
  par(mfrow = c(1,1))
  
  barchart(Percent_Ventilated ~ Day, 
           groups=Classifier, 
           ventData2, 
           auto.key = list(columns = 3),
           par.settings = list(
             superpose.polygon=list(col=c("#A3C1AD","darkblue"), border="black"),
             strip.background=list(col="red"),
             strip.border=list(col="black")),
           ylim = c(0,1),
           main="Percentage ventilated for each classifier")
  
  dev.off()
}
