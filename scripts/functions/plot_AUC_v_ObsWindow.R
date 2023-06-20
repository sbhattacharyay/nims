#### Supplementary function: Plot AUC vs. observation window ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk

plot.AUC.ObsWindow <- function(compiled.AUC.df,
                               ow.cutoff = 60, 
                               ow.units = 'min',
                               num.col=3,
                               axis.title.font.size = 8,
                               panel.title.font.size = 7,
                               axis.text.font.size = 5,
                               auc.max = 0.8,
                               auc.min = 0.2,
                               step.size = 3){
  
  if (ow.units == 'min' | ow.units == 'm' | ow.units == 'mins' ) {
    AUC.curves.plot <- compiled.AUC.df %>%
      mutate(ObsWindow = ObsWindow*60) %>%
      filter(ObsWindow <= ow.cutoff) %>%
      ggplot(aes(x = ObsWindow,y = meanValue)) + 
      facet_wrap(~Threshold,
                 scales = 'free',
                 ncol = num.col)+
      xlab("Observation Window (min)") + 
      ylab("AUC") +
      coord_cartesian(ylim = c(auc.min,auc.max)) + 
      geom_segment(x = 3, y = .5, xend = ow.cutoff, yend = 0.5,linetype = "dashed",size=.75/.pt, color = 'black')+
      geom_ribbon(aes(ymin=lowerValues,ymax=upperValues),alpha = 0.1,fill='blue') +
      geom_line(color = 'blue',size=1.3/.pt) + 
      geom_pointrange(aes(ymin=lowerValues,ymax=upperValues),fatten=1,size=1/.pt) +
      scale_y_continuous(breaks = seq(auc.min,auc.max,by = 0.1)) + 
      scale_x_continuous(breaks = seq(3,ow.cutoff,by = step.size)) + 
      theme_bw()+
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=panel.title.font.size, color = "black",face = 'bold'), 
        axis.text.x = element_text(size = axis.text.font.size, color = "black"),
        axis.text.y = element_text(size = axis.text.font.size, color = "black"),
        axis.title.x = element_text(size = axis.title.font.size, color = "black",face = 'bold'),
        axis.title.y = element_text(size = axis.title.font.size, color = "black",face = 'bold'),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
      )
  } else {
    AUC.curves.plot <- compiled.AUC.df %>%
      filter(ObsWindow <= ow.cutoff & (ObsWindow == 0.05 | ObsWindow >= 0.5)) %>%
      ggplot(aes(x = ObsWindow,y = meanValue)) + 
      facet_wrap(~Threshold,
                 scales = 'free',
                 ncol = num.col)+
      xlab("Observation Window (hr)") + 
      ylab("AUC") +
      coord_cartesian(ylim = c(auc.min,auc.max)) + 
      geom_segment(x = 0.05, y = .5, xend = ow.cutoff, yend = 0.5,linetype = "dashed",size=.75/.pt, color = 'black')+
      geom_ribbon(aes(ymin=lowerValues,ymax=upperValues),alpha = 0.1,fill='blue') +
      geom_line(color = 'blue',size=1.3/.pt) + 
      geom_pointrange(aes(ymin=lowerValues,ymax=upperValues),fatten=1,size=1/.pt) +
      scale_y_continuous(breaks = seq(auc.min,auc.max,by = 0.1)) + 
      scale_x_continuous(breaks = seq(0,ow.cutoff,by = step.size)) + 
      theme_bw()+
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=panel.title.font.size, color = "black",face = 'bold'), 
        axis.text.x = element_text(size = axis.text.font.size, color = "black"),
        axis.text.y = element_text(size = axis.text.font.size, color = "black"),
        axis.title.x = element_text(size = axis.title.font.size, color = "black",face = 'bold'),
        axis.title.y = element_text(size = axis.title.font.size, color = "black",face = 'bold'),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
      )
  }
  return(AUC.curves.plot)
}