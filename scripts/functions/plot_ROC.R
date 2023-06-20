#### Supplementary function: Plot receiver operating characteristic (ROC) curves ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk

plot.ROC <- function(plot.roc.axes.df,
                     num.col=3,
                     axis.title.font.size = 8,
                     panel.title.font.size = 7,
                     axis.text.font.size = 5){
  
  ROC.curves.plot <- plot.roc.axes.df %>% 
    ggplot(aes(x = FPR)) +
    facet_wrap( ~ Threshold,
                scales = 'free',
                ncol = num.col) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    coord_cartesian(ylim = c(0,1),xlim = c(0,1))+
    geom_segment(x = 0, y = 0, xend = 1, yend = 1,alpha = 0.5,linetype = "dashed",size=.75/.pt, color = 'gray')+
    geom_ribbon(aes(ymin = lowerTPR, ymax = upperTPR), alpha = 0.1,fill='red',linetype = "dotdash",size=.75/.pt,color='black') +
    geom_line(aes(y = meanTPR), alpha = 1, size=1.3/.pt,color='red') +
    theme_classic()+
    theme(
      strip.text = element_text(size=panel.title.font.size, color = "black",face = 'bold'), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = axis.text.font.size, color = "black"),
      axis.text.y = element_text(size = axis.text.font.size, color = "black"),
      axis.title.x = element_text(size = axis.title.font.size, color = "black",face = 'bold'),
      axis.title.y = element_text(size = axis.title.font.size, color = "black",face = 'bold'),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
      aspect.ratio = 1,
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
}
  