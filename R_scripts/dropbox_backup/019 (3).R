disc <- c(5, 3, 0, 2, 0, 3, 2, 3, 6, 1, 2, 1, 2, 1, 3, 3, 3, 5, 2, 4,
4, 0, 2, 3, 7, 12, 3, 10, 9, 2, 3, 7, 7, 2, 3, 3, 6, 2, 4, 3, 5, 2, 
2, 4, 0, 4, 2, 5, 2, 3, 3, 6, 5, 8, 3, 6, 6, 0, 5, 2, 2, 2, 6, 3, 4, 
4, 2, 2, 4, 7, 5, 3, 3, 0, 2, 2, 2, 1, 3, 4, 2, 2, 1, 1, 1, 2, 1, 4, 
4, 3, 2, 1, 4, 1, 1, 1, 0, 0, 2, 0);
tbl <- table(disc)
nm <- names(tbl)

# funkcija "table" prebroji koliko se pojedini element (kategorija)
# pojavljuje u vektoru i napravi tablicu u kojoj je svakom
# jedinstvenom elementu pridružen broj ponavljanja u vektoru. 

# tbl varijabla sadrži tablicu koja pokazuje koliko je puta u 
# godinama od 1850. do 1959. zabilježen odreðen broj znaèajnih 
# otkriæa. Npr. u tom je razdoblju u 12 pojedinaènih godina
# zabilježeno 1 otkriæe, u 26 pojedinèanih godina 2 otkriæa itd.  

# nm varijabla sadrži "imena" podataka u tablici, u našem sluèaju
# imena su brojevi otkriæa (0 otkiæa, 1 otkriæe, 2 otkriæa itd.).
# nm je vektor koji sadrži character data type pa se podaci u njemu
# ne mogu manipulirati brojèanim operatorima (+ - * /). Za dijeljenje
# svakog èlana nm vektora sa dva podatke u njemu treba prije pretvoriti
# u numerièki tip podataka funkcijom as.numeric(nm):

nm_num <- as.numeric(nm)
nm_num <- nm_num/2
nm_num

# raèunanje prosjeènog broja bitnih otkiæa u periodu 1860-1959:
tbl_numeric <- as.numeric(tbl)
tbl_numeric
mean(tbl_numeric)
