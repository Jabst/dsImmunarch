geneUsageDS <- function (dataframe, gene) {
  library(dsBase)
  
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  if(!isValidDS(dataframe)) {
    return (NA)
  }
  
  library(immunarch)
  ret <- parse_mixcr(dataframe)
  return(geneUsage(ret, .gene = gene))
}