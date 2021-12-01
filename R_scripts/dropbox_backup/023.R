lst <- list(vc1 = c(1, 3, 5), vc2 = c(2, 4, 6), vc3 = c(10, 20, 30, 40, 50))
lst

lst$vc4 <- c("a", "b", "xyz")
lst
  # sada je na listi lst vektor vc4

m1 <- matrix(1:30, nrow = 5, ncol = 6, byrow = TRUE)
m1

m2 <- matrix(1:6, nrow = 2, ncol = 3)
m2

m1 [4:5, 4:6] <- m2
m1

lst$m1 <-  m1
lst
  # sada je na listi lst matrica m1

b <- lst["vc3"]
b
d <- lst[["vc3"]]
d
class(b)
class(d)

  # Kada pozovemo vc3 pomoæu lst[["vc3"]], vraæa nam se vrijednost vc3 s liste
  # kao vektor brojeva ( class(lst[["vc3"]]) je "numeric"). 
  # Kada vc3 pozovemo pomoæu lst["vc3"], vraæa nam se lista duljine 1 koja 
  # sadri samo vektor vc3, class(lst["vc3"]) je "list". 
  # "[" vraæa listu, a "[[" vraæa tip podataka koji je spremljen na listu
  # (vektor, matrix...)

a <- lst[1:3]
a
b <- lst[[1:3]]
  # pozivanje više elemenata liste je moguæe samo sa "[", a kad se koristi 
  # "[[" dolazi do errora.

e <- lst$vc3
e
