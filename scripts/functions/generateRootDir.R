# Generates standardised directory for graph plotting which updates based on the date

generateRootDir <- function() {
  rootDir <- file.path("../plots",Sys.Date())
  dir.create(rootDir, showWarnings = FALSE, recursive = TRUE)
  return(rootDir)
}