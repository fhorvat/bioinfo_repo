m1 <- matrix(1:30, nrow = 5, ncol = 6, byrow = TRUE)
m1

m2 <- matrix(1:6, nrow = 2, ncol = 3)
m2

m1 [4:5, 4:6] <- m2
m1

mx1 <- m1
mx1 <- mx1[, -2]
mx1
