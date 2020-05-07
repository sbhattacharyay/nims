get_precrec_info <- function(list_of_preds){
  
  outList <- (lapply(seq_along(list_of_preds), function(x) {
    lapply(seq_along(list_of_preds[[1]][[1]])[-1], function(y) {
      s<-join_scores(list_of_preds[[x]][[1]][y],
                     list_of_preds[[x]][[2]][y],
                     list_of_preds[[x]][[3]][y],
                     list_of_preds[[x]][[4]][y],
                     list_of_preds[[x]][[5]][y], chklen = FALSE)
      l<-join_labels(list_of_preds[[x]][[1]][1],
                     list_of_preds[[x]][[2]][1],
                     list_of_preds[[x]][[3]][1],
                     list_of_preds[[x]][[4]][1],
                     list_of_preds[[x]][[5]][1], chklen = FALSE)            
      evalmod(scores=s,labels=l, dsids=1:5)})
  }))
  names(outList) <- names(list_of_preds)
  lapply(seq_along(outList), function(x){names(outList[[x]]) <<- names(list_of_preds[[x]][[1]])[-1]})
  return(outList)
}