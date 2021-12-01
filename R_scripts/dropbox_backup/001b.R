kocka_PMF <- function (side){
  d <- matrix (nrow=6, ncol=2)
  d[,1] <- 1:6
  d[,2] <- 1/6
  dimnames(d) <- list (c("a","b","c","d","e","f"), c("ishod", "vjerojatnost"))
  if (side %in% d[,1]){
    return (1/6)
  }
  else {
    warning ("Argument is non-discrete")
    return (0)
  }
}
kocka_PMF <- Vectorize (kocka_PMF)
kocka_PMF (c(12,2,3,4,0.6,0.2,78))


