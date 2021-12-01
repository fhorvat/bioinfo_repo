rm(list=ls())
# my.rnorm
my.rnorm <- function (n, mean, sd){
  set.seed(2)
  b <- runif(n)
  d <- qnorm(b, mean, sd)
  return(d)
}
my.rnorm(3, 0, 1)

# test
set.seed(2)
rnorm(3, 0, 1)
# my.rnorm vraæa neke iste brojeve kao i rnorm za isti n 
# (uz isti seed), ali neki brojevi nisu isti, ne znam zašto :)


rm(list=ls())
# my.rpois
my.rpois <- function (n, l){
  set.seed(3)
  b <- runif(n)
  d <- qpois(b, l)
  return(d)
}
my.rpois(13, 3)

# test
set.seed(3)
rpois(13, 3)
# my.rpois vraæa iste brojeve kao i rpois za isti n i lamda 
# (uz isti seed)


rm(list=ls())
# my.runif
my.runif <- function (n, min, max){
  a <- seq(min, max, length = 20000000)
  set.seed(3)
  b <- sample(a, n)
  return(b)
}
my.runif(3, 1, 5)

# test
set.seed(3)
runif(3, 1, 5)
# my.runif vraæa iste brojeve kao i runif za isti n, min i max 
# (uz isti seed)
