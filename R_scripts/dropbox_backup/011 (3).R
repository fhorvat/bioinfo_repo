sinus <- function (x=c(0,2*pi)){
  x <- seq(x[1], x[2], 0.01)
  plot(x, sin(x), type="l", xlim=c(x[1], x[2]))
}
sinus(c(1,10))


