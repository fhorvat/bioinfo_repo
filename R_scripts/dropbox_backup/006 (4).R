manh <- function(v1,v2){
  return(sum(abs(v1-v2)))
}

eucl <- function (v1,v2){
  return(sqrt(sum((v2-v1)^2)))
}

cheb <- function (v1,v2){
  return(max(abs(v2-v1)))
}

dist <- function (v1,v2,metric){
  if (metric == 'manhattan'){manh (v1,v2)} else 
      if (metric == 'euclidean') {eucl (v1,v2)} else 
        if (metric == 'chebyshev') {cheb (v1,v2)} else 
        print ('Error: distance metric must be manhattan, euclidean or chebyshev')
}
metric <- ('ManHattan')
v1 <- c(4,9,18)
v2 <- c(5,7,4)
metric <- tolower(metric)
dist (v1,v2,metric)

