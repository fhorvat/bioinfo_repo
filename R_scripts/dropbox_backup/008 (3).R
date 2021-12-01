kut <- function (v1,v2){
  acos(sum(v1*v2)/(sqrt(sum(v1*v1))*sqrt(sum(v2*v2))))
}
v1 <- c(3,0)
v2 <- c(5,5)
rad <- kut (v1,v2)
deg <- rad*(180/pi)
deg

 
