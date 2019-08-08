repOverlapDS <- function (dataframe) {
  library(dsBase)
  
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  check <- isValidDS(dataframe)
  
  # return missing value if the input vector is not valid
  if (!check) {
    return (NA)
  }
  
  library(immunarch)
  ret <- parse_mixcr(dataframe)
  
  overlapResult = repOverlap(ret, .method = "morisita", .verbose = F)
  
  return (overlapResult)
}