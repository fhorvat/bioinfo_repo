class(iris)
iris

# class objekta iris je data.frame


a <- nrow(iris)
a
b <- ncol(iris)
b

# broj redova data.frame-a se može saznati funkcijom nrow, a broj 
# stupaca faunkcijom ncol. Broj redova je 150, a broj stupaca 5. 


sapply(iris, class)

# class varijabli unutar data.frame-a možemo saznati tako da
# funkciju class primjenimo na cijeli data.frame pomoæu funkcije
# sappply. Rezultat je popis klasa varijabli u svakom stupcu.  


petal.width1 <- iris[[10,4]]
petal.width1

petal.width2 <- iris$Petal.Width[10]
petal.width2

# petal.width1 je vrijednost širine latice u 10. redu dobivena 
# pomoæu "[[", a petal.width2 je vrijednost širine latice u 10. redu
# dobivena pomoæu "$"


petal.width3 <- iris$Petal.Width
f.petal.width3 <- factor(petal.width3)
l.f.petal.width3 <- levels (f.petal.width3)
l.f.petal.width3

# l.f.petal.width2 pokazuje sve jedinstvene vrijednosti širine
# latica, a dobiven je tako da je prvo iz data.frame iris u vektor 
# izdvojen cijeli stupac koji pokazuje širinu latica (petal.width3).
# Zatim je taj vektor pretvoren u faktor (f.petal.width3) i na kraju
# su iz faktora izdvojeni leveli vrijednosti tog faktora 
# (l.f.petal.width3) koji prikazuju jedinstvene vrijednosti podataka
# unutar faktora. 


virginica <- iris[iris$Species == "virginica", ]
virginica

# virginica je data.set koji sadrži podatke iz iris samo za 
# vrstu I. virginica


sepal.width.versicolor1 <- iris[iris$Species == "versicolor", ]
sepal.width.versicolor2 <- sepal.width.versicolor1$Sepal.Width
sepal.width.versicolor2

# širina lapova vrste I. versicolor nalaze se u vektoru 
# sepal.width.versicolor2 
 

my.iris <- iris
my.iris.petal.width.sqrd <- sqrt(my.iris$Petal.Width) 
my.iris$Petal.Width.Sqrd <- my.iris.petal.width.sqrd
my.iris

# my.iris sada ima 6. stupac Petal.Width.Sqrd u kojem su vrijednosti
# korijena širine latica (stupac Petal.Width)


my.iris$Sepal.Length <- NULL
my.iris

# my.iris sada više nema stupac Sepal.Length
