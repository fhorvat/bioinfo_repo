f <- function (t,x){
  n <- 0:t
  zbr <- (((-1)^n)/(factorial((2*n)+1)))*(x^((2*n)+1))
  y <- sum(zbr)
    return(y)
}

x <- seq(from=-5, to=5, length.out = 100)
y1 <- 0;
y2 <- 0;
y3 <- 0;
for (i in seq_along(x) ){
  y1[i] <- f(1, x[i]);
  y2[i] <- f(3, x[i]);
  y3[i] <- f(5, x[i]);
}
plot (x, y1, type="l")
lines (x, y2)
lines (x, y3)

