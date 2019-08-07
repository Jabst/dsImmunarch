kmers <- function (dataframe) {
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  check <- isValidDS(dataframe)
  
  library(immunarch)
  ret <- parse_mixcr(dataframe)
  
  # return missing value if the input vector is not valid
  #if(!check) {
    kmers = getKmers(dataframe, 10, .col = "aa")
    result <- kmer_profile(kmers, "prob")
  #} else {
    result <- NA
  #}
  
  return(result)
}