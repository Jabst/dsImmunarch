repClonality <- function (dataframe) {
  library(dsBase)
  
  print(class(dataframe))
  
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  check <- isValidDS(dataframe)
  
  library(immunarch)
  ret <- parse_mixcr(dataframe)
  
  # return missing value if the input vector is not valid
  if(!check) {
    result <- repClonality(ret, .gene = "clonal.prop")
  } else {
    result <- NA
  }
  
  return(result)
}