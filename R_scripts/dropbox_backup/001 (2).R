# infix operatori su funkcije koje se pozivaju tako da ime funkcije dolazi 
# izmeðu dva argumenta (za razliku od veæine funkcija koje su "prefix" operatori
# i kod kojih ime funkcije dolazi ispred argumenata). 

# %in% provjerava postoji li lijevi argument u desnom i vraæa
# logièku vrijednost (T ako postoji, F ako ne postoji). Primjer:

x <- 1:10
2 %in% x # vraæa TRUE jer 2 postoji u x
11 %in% x # vraæa FALSE jer 11 ne postoji u x
"b" %in% x # vraæa FALSE je "b" ne postoji u x

