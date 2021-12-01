x=seq(2,12,length=200)
PDF <- function(x){
  dunif (x, min=2, max=12)
}
CDF <- function (x){
  punif(x, min=2, max=12)
}

plot(x,PDF(x),type="l",xlim=c(1,12),ylim=c(0,1),
     lwd=1,col="red")
lines(x,CDF(x),type="l",lwd=1,col="blue")