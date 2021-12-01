par(mfrow=c(2,1), mar=c(2,2,2,1))
#plot1
mean1 <- 0
sd1 <- 2
a <- mean1+sd1
b <- mean1-sd1
x <- seq(-4,4,length=200)*sd1+mean1
y <- dnorm(x,mean1,sd1)
plot(x, y, type="l", bty="l")
i <- x >= b & x <= a
polx <- c(b,x[i],a)
poly <- c(0,y[i],0)
polygon(polx,poly, col="red")

#plot 2
mean2 <- 3
sd2 <- 1
a2 <- mean2+sd2
b2 <- mean2-sd2
x2 <- seq(-4,4,length=200)*sd2+mean2
y2 <- dnorm(x2,mean2,sd2)
plot(x2, y2, type="l", bty="l")
i2 <- x2 >= b2 & x2 <= a2
polx2 <- c(b2,x2[i2],a2)
poly2 <- c(0,y2[i2],0)
polygon(polx2,poly2, col="green")