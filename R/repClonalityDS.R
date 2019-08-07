repClonalityDS <- function (dataframe) {
  library(dsBase)
  
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  check <- isValidDS(dataframe)
  
  # return missing value if the input vector is not valid
  if (!check) {
    library(immunarch)
    ret <- parse_mixcr(dataframe)
    return (immunarch::repClonality(ret, .gene = "clonal.prop"))
  } else {
    return (NA)
  }
}