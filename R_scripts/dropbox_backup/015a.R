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
  # (uz isti seed), a neki brojevi nisu isti, ne znam zašto :)

