multinomial_to_ordinal <- function(labels){
  no.classes <- length(levels(labels))
  
  ordinal.labels <- as.data.frame(matrix(nrow = length(labels), ncol = no.classes-1))
  
  for (m in 1:ncol(ordinal.labels)){
    ordinal.labels[,m] <- as.factor(as.numeric(as.numeric(labels) > m))
  } 
  
  return(ordinal.labels)
}