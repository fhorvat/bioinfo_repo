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
