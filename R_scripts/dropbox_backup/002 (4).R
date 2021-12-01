test.integer <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if (test==TRUE){ 
    return(TRUE)
  }
  else{ 
    return(FALSE) 
  }
}

kocka_CDF<-function(x){
  if (test.integer(x)==FALSE){
    warning ("Argument is non-discrete")
    return (0)
  }
  if (x>=1 & x<=6){
    return(floor(x)/6)
  }
  else {
    return(0)
  }
}

kocka_CDF <- Vectorize (kocka_CDF)
kocka_CDF (c(1,2,4))