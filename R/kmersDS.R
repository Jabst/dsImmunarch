kmersDS <- function (dataframe) {
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  check <- isValidDS(dataframe)
  
  print(1)
  
  # return missing value if the input vector is not valid
  if (!check) {
    return (NA)
  }
  
  print(2)
  
  library(immunarch)
  ret <- parse_mixcr(dataframe)
  
  print(3)
  
  kmers = getKmers(dataframe, 10, .col = "aa")
  print(4)
  kmersProfie = kmer_profile(kmers, "prob")
  
  print(5)
  
  return (kmersProfie)
}