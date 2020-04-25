cross_val_splits <- function(dataTable,k){
  set.seed(123)
  ratio<-nrow(dataTable)/k
  permutation<-sample(1:nrow(dataTable))
  
  if (nrow(dataTable)%%k == 0){
    return(split(permutation, ceiling(seq_along(permutation)/ratio)))
  }else{
    idxV <- 1:((k-1)*(floor(ratio)))
    tempList<-split(permutation[idxV], ceiling(seq_along(permutation[idxV])/floor(ratio)))
    tempList[[as.character(k)]]<-permutation[(tail(idxV,1)+1):tail(nrow(dataTable),1)]
    return(tempList)
  }
}