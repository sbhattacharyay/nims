get_auc_plots <- function(auc_info){
  val<-names(auc_info[[1]])
  col_vector=c('red','green','blue')
  name_vector=c('CM: ','CF only: ', 'MF only: ')
  plot_sep<-list(c(1,4,6),c(2,4,7),c(3,5,8))
  plots <- vector(mode = "list", length = length(val))
  i = 1
  for (z in plot_sep){
    auc_stats<<-sapply(seq_along(z), function(x) {
                    sapply(seq_along(auc_info[z][[1]]), function(y) {
                      paste(name_vector[x],round(mean(auc_info[z][[x]][[y]]$fold.AUC),2), "(+/-",round(sd(auc_info[z][[x]][[y]]$fold.AUC),2),")")
                    })
    })
    lapply(seq_along(auc_info[z][[1]]), function(y){
      plots[[y]][[i]]<<-arrangeGrob((as.ggplot(function()
        sapply(seq_along(auc_info[z]), function(x) {
          plot(auc_info[z][[x]][[y]]$perf,col='grey',lty=3)
          plot(auc_info[z][[x]][[y]]$perf,avg="vertical",col=col_vector[x], lwd=2, add=TRUE)
          par(new=T)
        })
      )+coord_fixed()),top=textGrob(val[y],vjust = 1.5, gp = gpar(fontface = "bold", cex = 1.5)),
        bottom=textGrob(paste('Mean AUC (+/-1 SD)',paste(auc_stats[y,], collapse='\n'), sep='\n')))
      par(new=F) 
    })
    i=i+1
  }
  return(plots)
}