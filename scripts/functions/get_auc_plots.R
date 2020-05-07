get_auc_plots <- function(auc_info){
  val<-names(auc_info[[1]])
  col_vector=c('red','darkgreen','blue')
  plot_sep<-list(c(1:3),c(4:5),c(6:8))
  plots <- vector(mode = "list", length = length(val))
  i = 1
  for (z in plot_sep){
    lapply(seq_along(auc_info[z][[1]]), function(y){
      plots[[y]][[i]]<<-arrangeGrob(as.grob(function()
        sapply(seq_along(auc_info[z]), function(x) {
          plot(auc_info[z][[x]][[y]]$perf,col='grey',lty=3)
          if (i==1&y==1|i==3&y==1) legend(0.6,0.4, c("GOSE", "Fav","Death"),fill=col_vector, bty = "n",pt.cex = 2, y.intersp=2.2)
          if (i==2&y==1) legend(0.6,0.4, c("GOSE","Death"),fill=col_vector[-2], bty = "n",pt.cex = 2, y.intersp=2.2)  
          if (i==2) plot(auc_info[z][[x]][[y]]$perf,avg="vertical",col=col_vector[-2][x], lwd=2, add=TRUE)
          else plot(auc_info[z][[x]][[y]]$perf,avg="vertical",col=col_vector[x], lwd=2, add=TRUE)
          par(new=T)
        })
      ),top=textGrob(val[y],vjust = 1.5, gp = gpar(fontface = "bold", cex = 1.5)))
      par(new=F) 
    })
    i=i+1
    }
  return(plots)
}