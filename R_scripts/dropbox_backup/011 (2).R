a <- 0
b <- 0
d <- 0
a <- ppois(3, lambda = 6, lower=TRUE) # vjerojatnost za 3 ili manje
b <- ppois(4, lambda = 6, lower=FALSE) # vjerojatnost za 5 ili više
d <- 1-(b+a)
d # vjerojatnost za toèno 4

a <- 0
a <- ppois(8, lambda = 6, lower=FALSE) 
a # vjerojatnost za 9 ili više

x <- 0:25
y1 <- ppois(x, lambda = 6, lower=FALSE)
y2 <- ppois(x, lambda = 9, lower=FALSE)
y3 <- ppois(x, lambda = 12, lower=FALSE)
plot(x, y1, type = "l", col = "red", bty = "l")
lines(x, y2, col = "green")
lines(x, y3, col = "blue")

i <- x > 8
px <- x[i]
py1 <- y1[i]
py2 <- y2[i]
py3 <- y3[i]
points(px, py1, col="red", pch=2)
points(px, py2, col="green", pch=7)
points(px, py3, col="blue", pch=8)

# crvene, zelene i plave toèke oznaèavaju regiju na grafa koja 
# pokazuje vjerojatnost da se dionicama trgovalo 9 ili više puta 
# u minuti. Crvene su za dionice kojima se prosjeèno trguje 6 puta
# u minuti, zelene za dionice kojima se prosjeèno trguje 9 puta u 
# minuti, a plave su za dionice kojima se prosjeèno trguje 12 puta u
# minuti




