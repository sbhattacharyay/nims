get_precrec_plots <- function(precrec_info){
  col_vector=c('GOSE','Fav','Death')
  val<-names(precrec_info[[1]])
  plots <<- vector(mode = "list", length = length(val))
  plot_sep<-list(c(1:3),c(4:5),c(6:8))
  n = 1
  for (z in plot_sep){
    for (i in seq_along(precrec_info[z][[1]])){
    p_roc <<- ggplot()
    for (j in seq_along(precrec_info[z])){
      p_roc <<- p_roc + geom_line(data=cbind(Predictions=col_vector[j],subset(fortify(precrec_info[z][[j]][[i]]), curvetype == "PRC")),
                                  mapping=aes(x = x, y = y, color = Predictions))
    }
    if (i==1 & n==1 | i==1 & n==3) {
      p_roc <<- p_roc + xlab("Recall") + ylab("Precision") + scale_color_manual(name="", 
                         labels = c('GOSE','Fav','Death'), 
                         values = c('red', 'darkgreen', 'blue'))+
                theme(text=element_text(size=16),axis.text.x=element_text(hjust=1),panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),panel.background = element_blank(), 
                      axis.line = element_line(colour = "black"),legend.text=element_text(size=10))
    }
    
    else if (i==1 & n==2) {
      p_roc <<- p_roc + xlab("Recall") + ylab("Precision") + scale_color_manual(name="", 
                                                                                labels = c('GOSE','Death'), 
                                                                                values = c('red', 'blue'))+
        theme(text=element_text(size=16),axis.text.x=element_text(hjust=1),panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),legend.text=element_text(size=10))      
    }
    else if (i!=1 & n==2) {
      p_roc <<- p_roc + xlab("Recall") + ylab("Precision") + scale_color_manual(name="", 
                                                                                labels = c('GOSE','Death'), 
                                                                                values = c('red', 'blue'))+
        theme(legend.position = "none",text=element_text(size=16),axis.text.x=element_text(hjust=1),panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),panel.background = element_blank(), 
              axis.line = element_line(colour = "black"))      
    }
    else{
      p_roc <<- p_roc + xlab("Recall") + ylab("Precision")+scale_color_manual(name="", 
                             labels = c('GOSE','Fav','Death'), 
                             values = c('red', 'darkgreen', 'blue'))+
          theme(legend.position = "none",text=element_text(size=16),axis.text.x=element_text(hjust=1),panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),strip.background = element_blank(), 
                axis.line = element_line(colour = "black"))
    }

    plots[[i]][[n]] <<- arrangeGrob(p_roc,top=textGrob(val[i],vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
    }
    n=n+1
  }
  return (plots)
}