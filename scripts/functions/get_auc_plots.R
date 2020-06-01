get_auc_plots <- function(auc_info, auc_ci_info){
  val<-names(auc_info[[1]])
  col_vector=c('red','green','blue')
  name_vector=c('CM: ','CF ONLY: ', 'MF ONLY: ')
  plot_sep<-list(c(1,2,3))
  plots <- vector(mode = "list", length = length(val))
  i = 1
  for (z in plot_sep){
    # auc_stats<<-sapply(seq_along(z), function(x) {
    #                 sapply(seq_along(auc_info[z][[1]]), function(y) {
    #                   paste(name_vector[x],round(mean(auc_info[z][[x]][[y]]$fold.AUC),2), "(+/-",round(sd(auc_info[z][[x]][[y]]$fold.AUC),2),")")
    #                 })
    # })
    auc_stats<<-sapply(seq_along(z), function(x) {
      sapply(seq_along(auc_info[z][[1]]), function(y) {
        paste(name_vector[x],round(auc_ci_info[z][[x]][[y]]$cvAUC,2), "(",round(auc_ci_info[z][[x]][[y]]$ci[1],2),"-",
                                                                          round(auc_ci_info[z][[x]][[y]]$ci[2],2),")")
      })
    })    
    lapply(seq_along(auc_info[z][[1]]), function(y){
      plots[[y]][[i]]<<-arrangeGrob((as.ggplot(function()
        sapply(seq_along(auc_info[z]), function(x) {
          plot(auc_info[z][[x]][[y]]$perf,col='grey',lty=3)
          plot(auc_info[z][[x]][[y]]$perf,avg="vertical",col=col_vector[x], lwd=2, add=TRUE,
               show.spread.at = seq(0, 1, 0.3),spread.estimate="stderror",spread.scale=2,plotCI.col=col_vector[x],plotCI.lwd=2)
          par(new=T)
        })
      )+coord_fixed()),top=textGrob(val[y],vjust = 1.5, gp = gpar(fontface = "bold", cex = 1.5)),
        bottom=textGrob(paste('MEAN AUC (LOWER-UPPER)',paste(auc_stats[y,], collapse='\n'), sep='\n'),just=c(0.3,0)))
      par(new=F) 
    })
    i=i+1
  }
  return(plots)
}