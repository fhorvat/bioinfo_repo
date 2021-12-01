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
    return (0)
  }
  if (x>=1 & x<=6){
    return(1/6)
  }
  else {
    return(0)
  }
}

kocka_CDF<-function(x){
  if (test.integer(x)==FALSE){
    return (0)
  }
  if (x>=1 & x<=6){
    return(floor(x)/6)
  }
  else {
    return(0)
  }
}

x <- c(0.02, 0.1, 1,4 ,5,6,4,7)
kocka_PMF <- Vectorize (kocka_PMF)
kocka_CDF <- Vectorize (kocka_CDF)
matplot(x, cbind(kocka_CDF(x),kocka_PMF(x)),type="h",col=c("red","yellow"),
        lty=c(1,1),xlim=c(-1,7), ylim=c(0,1), bty="l", lwd=c(5,2))
axis(1, at=c(1,3,5,7))
