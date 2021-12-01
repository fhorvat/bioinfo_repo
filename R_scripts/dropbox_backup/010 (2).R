par(mfrow=c(2,1), mar=c(2,2,2,1))

#plot 1
#CDF
mean <- 0
sd <- 1
x <- seq(-4, 4, length.out = 1000) * sd + mean
y1 <- pnorm(x, mean, sd)
plot(x, y1, type = "l", bty = "l", col="red", ylim=c(0,1))
#f(x)=x
lines(x, x, col = "blue")
#qnorm 
y2 <- qnorm(x, mean, sd)
lines(x, y2, col="green")

rm(list=ls())

#plot 2
#CDF
mean <- 2
sd <- 3
x <- seq(-4, 4, length.out = 1000) * sd + mean
y1 <- pnorm(x, mean, sd)
plot(x, y1, type = "l", bty = "l", col="red", ylim=c(0,1))
#f(x)=x
lines(x, x, col = "blue")
#qnorm 
y2 <- qnorm(x, mean, sd)
lines(x, y2, col="green")

#qnorm i pnorm su meðusobno inverzne funkcije
x <- 0.5
a <- qnorm (x)
b <- pnorm (a)
all.equal(x,b)

qnorm(pnorm(3))
pnorm(qnorm(0.4))

#svejedno je napišemo li kompoziciju qnorm(pnorm(x)) ili 
#pnorm(qnorm(x)), u oba sluèaja vraæa se argument s kojim se
#ulazi u kompoziciju 