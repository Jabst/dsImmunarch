
#'
#' @title Computes statistical mean of a vectores
#' @description Calculates the mean value.
#' @details if the length of input vector is less than the set filter
#' a missing value is returned.
#' @param xvect a vector
#' @return a numeric, the statistical mean
#' @author Gaye, Not A.
#' @export
#'

#library(immunarch)

geneUsageDS <- function (arg1) {
  
  
  
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  # check <- isValidDS(xvect)
  
  library(immunarch)
  
  print("olá amigos")
  
  
  # return missing value if the input vector is not valid
  #if(!check){
    #result <- geneUsage(xvect, .gene = gene)
    result <- arg1
  #}else{
   # result <- NA
  #}
  
  return(result)
}