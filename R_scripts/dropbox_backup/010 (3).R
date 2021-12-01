# 0 to 5, 0.01 increment, w/o seq
n <- 1:501
a [n] <- 0+((n-1)*0.01)
a

#-pi to pi, w/o seq
n <- 1:25
a [n] <- -pi+((n-1)*((2*pi)/24))
a

#0 to 5, 0.01 increment, seq
a <- seq (0,5,0.01)
a

#-pi to pi, seq
a <- seq (-pi,pi,length.out=25)
a