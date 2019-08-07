repClonality <- function (dataframe) {
  library(dsBase)
  
  print(dim(dataframe))
  
  if(class(dataframe) == "character" | class(dataframe) == "integer" | class(dataframe) == "logical" | class(dataframe) == "numeric") {
    if(length(dataframe) > 0 & length(dataframe)  < 5) {
      
      print(1)
    } else {
      print(2)
    }
  }else{
    if(class(dataframe) == "factor"){
      tt <- tabulate(dataframe)
      xx <- which(tt > 0 & tt < 5)
      if(length(xx) > 0) {
        print(3)
      } else {
        print(4)
      }
    }else{
      if(class(dataframe) == "data.frame" | class(dataframe) == "matrix"){
        if(dim(dataframe)[1] > 0 & dim(dataframe)[1] < 5){
          print(5)
        }else{
          print(6)
        }
      }else{
        print(7)
      }
    }
  } 
  
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