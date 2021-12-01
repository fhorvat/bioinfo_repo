fact1 <- function (n) {
  if (n==0) return (1) else 
    return (n*fact1(n-1))
}
fact1(6)

fact2 <-function (n) {
  return (sqrt(2*pi*n)*((n/exp(1))^n))
}
fact2(4)
