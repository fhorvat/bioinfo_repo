# Ako u kvadarat imamo upisan krug polumjera R, tada je površina 
# kruga (R^2)*pi, a površina kvadrata je 4*(R^2). Ako podijelimo
# površinu kruga s površinom kvadrata dobijemo pi/4. Ako nasumièno
# "bacamo" toèke u taj kvadrat, odnos izmeðu površine kruga i 
# površine kvadrata možemo dobiti tako da "prebrojimo" toèke u krugu
# i izvan kruga. Pošto znamo iz geometrije da je taj odnos pi/4, 
# tada slijedi da je: 
# pi = 4*(broj toèaka u krugu/broj toèaka u kvadratu). Kako sve 
# toèke padaju u kvadrat (bilo unutar kruga, bilo izvan kruga),
# pi = 4*(broj toèaka u krugu/ukupan broj toèaka).

# U zadatku biramo brojeve od -1 do 1 po x i y osi, znaèi naš kvadrat
# je od -1 do 1 na x osi dugaèak i od -1 do 1 na y osi visok. 
# Krug koji mu je upisan ima polumjer 1. 

# Da bi vidjeli je li toèka unutar kruga ili izvan kruga koristimo
# Pitagorin pouèak. Uzmemo n toèaka s nasumiènim koordinatama x i y.
# Ako promatramo trokut kojemu su katete a i b udaljenost toèke 
# od x i y osi, tada  nepoznata stranica c tog trokuta odreðuje
# je li toèka unutar kruga ili izvan kruga. Ako je c <= 1 tada je
# toèka unutar kruga, a ako je c > 1 tada je toèka izvan kruga.  

My.pi <- function (n){
  a <- runif(n, min = -1, max = 1)
  b <- runif(n, min = -1, max = 1)
  c <- sqrt((a^2) + (b^2))
  inside <- sum(c <= 1)
  my.pi <- 4 * (inside / n)
  return (my.pi)
}
set.seed(1234)
My.pi(100000)

# plot
My.pi.plot <- function(iter){
  pi.plot <- 0
  for (i in 1:iter){
    pi.plot[i] <- My.pi(i)
  }
  plot(pi.plot[1:i], type = "l", ylim = c(0, 5), 
       col = "blue", xlab = "Number of iterations")
  abline(pi, 0, col = "red", lwd = 0.25)
}
set.seed(1234)
My.pi.plot(500)


