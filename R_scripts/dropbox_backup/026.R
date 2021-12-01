# diag zamjenjuje dijagonalu matrice, a može se koristiti i za dobivanje
# vrijednosti dijagonale matrice (diag(m), m je matrix)

mm1 <- matrix (1:9, nrow = 3, ncol = 3, byrow = TRUE)
mm1

diag (mm1) <- c(10, 20, 30)
mm1

# syntax je èudan, kao da je diag varijabla unutar matrice kojoj
# pridružujemo vrijednost.
# Google kaže: 
# diag(object) <- value je skraæeno od object <- "diag<-"(object, value)
# tj. diag (mm1) <- c(10, 20, 30) je"mm1 <- "diag<-"(mm1, c(10, 20, 30)). 

mm1 <- matrix (1:9, nrow = 3, ncol = 3, byrow = TRUE)
mm1 <- "diag<-"(mm1, c(10, 20, 30))
mm1

# Sve replacements funkcije u R-u su napravljene tako da mi pozivamo skraæeni
# oblik, a R ustvari preimenuje naš skraæeni oblik u duži oblik i zove 
# odreðenu replacements funkciju. Opæenito:
# FUN(obj) <- value 
# obj <- "FUN<-"(obj, value)

rm(list=ls())
"My.replacement<-" <- function(mtrx, value){
  b <- ceiling(nrow(mtrx)/2)
  if (nrow(mtrx) %% 2 == 0){
    mtrx[c(b), ] <- value
    mtrx[c(b+1), ] <- value
  }
  else {
    mtrx[b, ] <- value
  }
  return(mtrx)
}

mtrx <- matrix (1:9, nrow = 8, ncol = 4)
mtrx
My.replacement (mtrx) <- 100:103
mtrx
