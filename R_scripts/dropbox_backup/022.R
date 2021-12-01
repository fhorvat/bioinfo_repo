n <- 3
m <- 4
p <- 7

mt1 <- matrix(1:100, nrow = p, ncol = m)
mt1

mt2 <- matrix(1:100, nrow = n, ncol = p)
mt2

mt3 <- mt2 %*% mt1
a <- dim(mt3)
a
a[1] == n & a[2] == m

# Ako su dimenzije matrice mt1 m i n, dimenzije matrice mt3 dobivene 
# množenjem mt1 %*% mt2 m i p, tada su dimenzije matrice
# mt2 n i p. 

n <- 4
m <- 5
p <- 7
mt1 <- matrix(1:100, nrow = m, ncol = n)
mt2 <- matrix(1:100, nrow = n, ncol = p )
mt3 <- mt1 %*% mt2
a <- dim(mt3)
a[1] == m & a[2] == p


# Generalna pravila kod množenja:
# - m1 %*% m2 = m3: 
#   - dimenzije m1 su m i n
#   - tada su dimenzije m2 n i p
#   - tada su dimenzije m3 m i p

# - m2 %*% m1 = m3
#   - dimenzije m2 n i p
#   - tada su m1 su p i m
#   - tada su dimenzije m3 n i p

