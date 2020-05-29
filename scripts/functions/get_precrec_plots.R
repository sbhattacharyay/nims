get_precrec_plots <- function(precrec_info){
  col_vector=c('red', 'darkgreen', 'blue')
  val<-names(precrec_info[[1]])
  plot_sep<-list(c(1,2,3))
  plots <<- vector(mode = "list", length = length(val))
  n = 1
  for (z in plot_sep){
    for (i in seq_along(precrec_info[z][[1]])){
    p_roc <<- ggplot()
    for (j in seq_along(precrec_info[z])){
      p_roc <<- p_roc + geom_line(data=cbind(Predictions=col_vector[j],subset(fortify(precrec_info[z][[j]][[i]]), curvetype == "PRC")),
                                  mapping=aes(x = x, y = y, color = Predictions))
    }
      p_roc <<- p_roc + xlab("Recall") + ylab("Precision") + xlim(0,1)+ylim(0,1)+ theme_bw() +
                theme(legend.position='none',text=element_text(size=12),panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),panel.background = element_blank(), strip.background = element_blank(),
                      axis.line = element_line(colour = "black"),legend.text=element_text(size=10))

    plots[[i]][[n]] <<- arrangeGrob(p_roc+coord_fixed(),top=textGrob(val[i],vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
    }
    n=n+1
  }
  return (plots)
}

# scale_color_manual(name="", 
#                  labels = c('Combined Model','Clinical Feature Only','Motion Feature Only'), 
#                  values = c('red', 'darkgreen', 'blue'))+