test.integer <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if (test==TRUE){ 
    return(TRUE)
  }
  else{ 
    return(FALSE) 
  }
}

kocka_PMF<-function(x){
  if (test.integer(x)==FALSE){
    warning ("Argument is non-discrete")
    return (0)
  }
  if (x>=1 & x<=6){
    return(1/6)
  }
  else {
    return(0)
  }
}

kocka_PMF <- Vectorize (kocka_PMF)
kocka_PMF (c(1,3,9,4,-6))



