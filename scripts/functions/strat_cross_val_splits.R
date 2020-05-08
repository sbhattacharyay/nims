strat_cross_val_splits <- function(Y,k){
  set.seed(123)
  
  levelsY <- levels(Y)
  
  class1<-table(Y)[1]
  class2<-table(Y)[2]
  
  ratio1<- class1/k
  ratio2<- class2/k
  
  idxs1 <- which(Y == levelsY[1])
  idxs2 <- which(Y == levelsY[2])

  permutation1<-sample(idxs1)
  permutation2<-sample(idxs2)
  
  if (class1%%k == 0){
    tempList1<-split(permutation1, ceiling(seq_along(permutation1)/ratio1))
  }else{
    idxV <- 1:(k*(floor(ratio1)))
    tempList1<-split(permutation1[idxV], ceiling(seq_along(permutation1[idxV])/floor(ratio1)))
    for (i in 1:(class1%%k)){
      tempList1[[as.character(i)]]<-c(tempList1[[as.character(i)]],permutation1[length(permutation1)+1-i])
    }
  }
  
  if (class2%%k == 0){
    tempList2<-split(permutation2, ceiling(seq_along(permutation2)/ratio2))
  }else{
    idxV <- 1:(k*(floor(ratio2)))
    tempList2<-split(permutation2[idxV], ceiling(seq_along(permutation2[idxV])/floor(ratio2)))
    for (i in 1:(class2%%k)){
      tempList2[[as.character(k+1-i)]]<-c(tempList2[[as.character(i)]],permutation2[length(permutation2)+1-i])
    }
  }
  
  return(mapply(c, tempList1, tempList2, SIMPLIFY=FALSE)%>%lapply(sample))
}