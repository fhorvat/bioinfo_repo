#koordinatni sustav
x_os <- seq(-2, 2, 0.333)
y_os <- seq(-3, 3, 0.5)
plot(xy.coords (x_os, y_os), main="Exercise 16", type="n", xlab="x values",
     ylab="y values", bty="l")
grid(nx=NULL, ny = NULL, col = "lightgray", lty = "dotted",
    lwd = par("lwd"))
axis(2, at=-2)

#linija1
x_lin <- seq(-2, -0.5, 0.5)
y_lin <- seq(-2.5, -1, 0.5)
lines(range(x_lin), range(y_lin), lwd=3)

#toèke1
x_pnt1 <- seq(-0.5, 0.5, 0.25)
y_pnt1 <- seq(-1, 1, 0.5)
points(xy.coords (x_pnt1,y_pnt1), col="red")

#toèke2 
x_pnt2 <- seq(0.5, 2, 0.5)
y_pnt2 <- seq(1, 2.5, 0.5)
points(xy.coords(x_pnt2, y_pnt2), cex=0.75, pch=4)

#strelica
arrows(-1, 1, -1.25,-1.25, length=0.15)

#text
text(-1, 1.5, labels="Line width 3", cex=0.75)

