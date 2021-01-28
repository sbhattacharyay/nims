movingAverage <- function(x, n=1, fill.type = NA) {
  zoo.filter <- zoo::rollmean(x, k = n,fill = fill.type)
  return(zoo.filter)
}