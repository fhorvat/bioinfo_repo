mat <- matrix(rnorm(30, 5, 2), ncol = 6, nrow = 5)
mat

sd.row <- apply(mat, 1, sd)
sd.row

sd.col <- apply(mat, 2, sd)
sd.col
