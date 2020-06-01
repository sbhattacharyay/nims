plotBeanplots <- function(variable, plotName, plotUnits, type, inputData, legend = TRUE, titled = TRUE, splitName = "death", splitLevel1 = 0, splitLevel2 = 1) {
  
  #pull out variable i into a long dataframe
  tempDataFrame = gather(dplyr::select(inputData, all_of(variable)))
  tempDataFrame$splitter <- inputData[[splitName]]
  
  #remove variable from names of days and reconstruct naming
  
  level.1.DF <- filter(tempDataFrame, splitter == splitLevel1)
  level.2.DF <- filter(tempDataFrame, splitter == splitLevel2)
  
  #set up beanplot
  if (type == "time") {
    xlimits <- c(0,6)
  } else {
    xlimits <- c(0,2)
  }
  
  ylimits <- c(min(tempDataFrame$value, na.rm = TRUE), max(tempDataFrame$value, na.rm = TRUE))
  
  #set plotting units
  plotUnits <- ifelse(!is.na(plotUnits), paste(" / ", plotUnits, sep = ""), "")
  
  #plot
  beanplot(value~splitter, level.1.DF, 
           #side = "first", 
           col=c("#A3C1AD","#A3C1AD"),
           xlim = c(0,5),
           ylim = ylimits,
           main = ifelse(titled == TRUE,
                         paste("Beanplot showing distribution of",
                               plotName,
                               "by vital status at discharge"), NA),
           ylab = paste(plotName, plotUnits, sep = ""),
           log = "",
           what = c(1,1,0,1), 
           las = 1,
           bw="nrd0")
  
  #plot deceased/level2
  beanplot(value~splitter, level.2.DF, 
           #side = "second", 
           col=c("darkblue","darkblue", "lightgrey"), 
           add = TRUE,
           what = c(1,1,0,1),
           bw="nrd0")
  
  #add arithmetic mean lines
  if (type == "time") {
    for (day in 1:5) {
      meanVal <- mean(level.1.DF$value[level.1.DF$key == paste("Day ",day,sep="")], na.rm = TRUE)
      lines(x = c(day, day - 0.5), 
            y = c(meanVal,meanVal), 
            col = "black", 
            lwd = 2)
      
      meanVal <- mean(level.2.DF$value[level.2.DF$key == paste("Day ",day,sep="")], na.rm = TRUE)
      lines(x = c(day, day + 0.5), 
            y = c(meanVal,meanVal), 
            col = "black", 
            lwd = 2)
    } 
  } else {
    meanVal <- mean(level.1.DF$value, na.rm = TRUE)
    lines(x = c(1, 1 - 0.5), 
          y = c(meanVal,meanVal), 
          col = "black", 
          lwd = 2)
    
    meanVal <- mean(level.2.DF$value, na.rm = TRUE)
    lines(x = c(1, 1 + 0.5), 
          y = c(meanVal,meanVal), 
          col = "black", 
          lwd = 2)
  }
  
  #add a legend
  if(legend == TRUE) legend('topright', fill=c('#A3C1AD','darkblue'), legend= c(splitLevel1, ifelse(splitLevel2 == "dead", "deceased", splitLevel2)), bty ="n")
}