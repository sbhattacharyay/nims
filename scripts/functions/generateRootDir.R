# Generates standardised directory for graph plotting which updates based on the date

generateRootDir <- function(drctry="") {
  rootDir <- file.path("../plots",Sys.Date(),drctry)
  dir.create(rootDir, showWarnings = FALSE, recursive = TRUE)
  return(rootDir)
}