cross_val_splits <- function(dataTable,k){
  set.seed(123)
  ratio<-nrow(dataTable)/k
  permutation<-sample(1:nrow(dataTable))
  
  if (nrow(dataTable)%%k == 0){
    return(split(permutation, ceiling(seq_along(permutation)/ratio)))
  }else{
    idxV <- 1:(k*(floor(ratio)))
    tempList<-split(permutation[idxV], ceiling(seq_along(permutation[idxV])/floor(ratio)))
    for (i in 1:(nrow(dataTable)%%k)){
      tempList[[as.character(i)]]<-c(tempList[[as.character(i)]],permutation[length(permutation)+1-i])
    }
    return(tempList)
  }
}