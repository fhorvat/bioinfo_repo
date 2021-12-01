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
  # my.rpois vraæa iste brojeve kao i rpois za isti n, lamda 
  # (uz isti seed)