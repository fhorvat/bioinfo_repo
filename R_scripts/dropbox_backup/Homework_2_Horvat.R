#############
# 1. excerise
#############
rm(list=ls())
test.integer <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if (test==TRUE){ 
    return(TRUE)
  }
  else{ 
    return(FALSE) 
  }
}

kocka_PMF<-function(x){
  if (test.integer(x)==FALSE){
    warning ("Argument is non-discrete")
    return (0)
  }
  if (x>=1 & x<=6){
    return(1/6)
  }
  else {
    return(0)
  }
}

kocka_PMF <- Vectorize (kocka_PMF)
kocka_PMF (c(1,3,9,4,-6))


#############
# 2. excerise
#############
rm(list=ls())
test.integer <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if (test==TRUE){ 
    return(TRUE)
  }
  else{ 
    return(FALSE) 
  }
}

kocka_CDF<-function(x){
  if (test.integer(x)==FALSE){
    warning ("Argument is non-discrete")
    return (0)
  }
  if (x>=1 & x<=6){
    return(floor(x)/6)
  }
  else {
    return(0)
  }
}

kocka_CDF <- Vectorize (kocka_CDF)
kocka_CDF (c(1,2,4,7,0.1))

#############
# 3. excerise
#############
rm(list=ls())
test.integer <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if (test==TRUE){ 
    return(TRUE)
  }
  else{ 
    return(FALSE) 
  }
}

kocka_PMF<-function(x){
  if (test.integer(x)==FALSE){
    return (0)
  }
  if (x>=1 & x<=6){
    return(1/6)
  }
  else {
    return(0)
  }
}

kocka_CDF<-function(x){
  if (test.integer(x)==FALSE){
    return (0)
  }
  if (x>=1 & x<=6){
    return(floor(x)/6)
  }
  else {
    return(0)
  }
}

x <- c(0.02, 0.1, 1,4 ,5,6,4,7)
kocka_PMF <- Vectorize (kocka_PMF)
kocka_CDF <- Vectorize (kocka_CDF)
matplot(x, cbind(kocka_CDF(x), kocka_PMF(x)), type="h", 
        col=c("red","yellow"), lty=c(1,1), xlim=c(-1,7), ylim=c(0,1), 
        bty="l", lwd=c(5,2))
axis(1, at=c(1,3,5,7))

#############
# 4. excerise
#############
rm(list=ls())
x=seq(2,12,length=200)
PDF <- function(x){
  dunif (x, min=2, max=12)
}
CDF <- function (x){
  punif(x, min=2, max=12)
}

plot(x,PDF(x),type="l",xlim=c(1,12),ylim=c(0,1),
     lwd=1,col="red")
lines(x,CDF(x),type="l",lwd=1,col="blue")

#############
# 5. excerise
#############
# U prvom zadatku raèunamo mass probability function
# diskretne uniformne distribucije. Diskretna distribucija znaèi
# da brojevi pri bacanju kocke mogu zauzeti toèno unaprijed odreğene 
# vrijednosti i zato su na grafu predstavljeni toèkom. Pri bacanju 
# kocke ishod ne moe zauzeti vrijednost npr. 2.5 ili 3.8 nego samo
# cijeli broj izmeğu 1 i 6. Uniformna distribucija znaèi da postoji 
# jednaka vjerojatnost da ishod bacanja kocke bude bilo koji cijeli 
# broj izmeğu 1 i 6. 

# Dunif se koristi za izraèun probability density function
# kod kontinuirane uniformne distribucije. Kod kontinuirane 
# uniformne distribucije svi intervali jednake duljine izmeğu
# minimalne i maximalne vrijednosti su jednako vjerojatni. Zato na 
# grafu to vidimo kao crtu, tj. interval. 

# Kada bi u prvom zadatku koristili dunif ne bi raèunali 
# probability mass function bacanja kocke kao sada (koliko je 
# vjerojatno da ishod bacanja kocke bude npr. 2) nego probability 
# density function (koliko je vjerojatno da kocka zauzme bilo koju
# vrijednost izmeğu 0 i 6). Bez obzira koristimo li kao argument u
# funkciji 2 ili 4.5 dobiveni ishod je  0.1666667

rm(list=ls())
x <- dunif (2, max=6)
y <- dunif (4.5, max=6)
if (x==y) print("ishod je isti bez obzira na argument") else print("nešto sam zeznuo")

#############
# 6. excerise
#############
rm(list=ls())
x <- sample (20, 10000, replace = TRUE)
x

#############
# 7. excerise
#############
# Funkcije Xnorm koriste se za operacije vezane uz
# normalnu distribuciju. 

# Normalna distribucija je kontinuirana distribucija 
# vjerojatnosti koja govori kolika je vjerojatnost da 
# neka primjeæena realna vrijednost padne izmeğu bilo 
# koje dvije realne granice. 

# dnorm je funkcija koja uzima argument x i vraæa visinu
# distribucije vjerojatnosti za svaku toèku unutar x. 
# Ako se koristi samo kao dnorm(x) uzima srednju vrijednost 0 
# i standardnu devijaciju 1:
rm(list=ls())
x <- seq(-4,4,length=200)
y <- dnorm (x)
plot(x,y,type="l",col="red")

# Moguæe je i definirati proizvoljnu srednju vrijednost i standardnu
# devijaciju: dnorm (x, mean=, sd=). Primjer:
rm(list=ls())
x <- seq(-10,10,length=200)
y <- dnorm (x, mean=2, sd=2)
plot(x,y,type="l",col="red", xlim=c(-10,10))

# Ukoliko standardna devijacija (sd) ostaje nepromjenjena, 
# poveæavanjem srednje vrijednosti (mean) dolazi do horizontalnog 
# pomicanja krivulje udesno, a smanjenje mean pomièe krivulju 
# horizontalno ulijevo:
rm(list=ls())
x <- seq(-15,15,length=200)
y1 <- dnorm (x, mean=0, sd=2)
y2 <- dnorm (x, mean=4, sd=2)
y3 <- dnorm (x, mean=-4, sd=2)
plot (x,y1,type="l",col="red", xlim=c(-15,15))
lines (x, y2, col="blue")
lines (x, y3, col="green")

# Smanjenje standardne devijacije (s) dovodi do uske i visoke 
# krivulje, a povišenje standardne devijacije ( s ) do široke 
# i niske krivulje: 
rm(list=ls())
x <- seq(-15,15,length=200)
y1 <- dnorm (x, mean=0, sd=4)
y2 <- dnorm (x, mean=0, sd=6)
y3 <- dnorm (x, mean=0, sd=2)
plot (x,y1,type="l",col="red", xlim=c(-15,15), 
      ylim=range(y1,y2,y3))
lines (x, y2, col="blue")
lines (x, y3, col="green")

# pnorm je funkcija kumulativne distribucije koja za vrijednost
# x daje vjerojatnost da je nasumièan broj dobiven normalnom 
# distribucijom manji od x. Prihvaæa iste argumente kao dnorm:
rm(list=ls())
x <- seq(-15,15,length=200)
y1 <- pnorm (x, mean=2, sd=4)
y2 <- pnorm (x)
plot (x,y1,type="l",col="red", xlim=c(-15,15), ylim=c(0,1))
lines(x,y2,col="green")

# Ako elimo da funkcija daje vjerojatnost da je nasumièan broj 
# dobiven normalnom distribucijom VEÆI od x moemo korisiti 
# argument lower.tail=FALSE:
rm(list=ls())
x <- seq(-15,15,length=200)
y1 <- pnorm (x, mean=2, sd=4, lower.tail=FALSE)
y2 <- pnorm (x, lower.tail=FALSE)
plot (x,y1,type="l",col="red", xlim=c(-15,15), ylim=c(0,1))
lines(x,y2,col="green")

# qnorm je funkcija inverzna funkciji pnorm: vrijednost argumenta
# x je vjerojatnost, a funkcija vraæa broj èija kumulativna 
# distribucija odgovara toj vjerojatnosti x. qnorm takoğer prihvaæa
# neobavezne argumente mean i sd:
rm(list=ls())
x <- seq(0,1,length.out=100)
y <- qnorm(x)
y2 <- qnorm(x, mean=2, sd=3)
plot(x,y, type="l",col="red", ylim=c(-10,10))
lines (x,y2, col='blue')

# rnorm je funkcija koja generira nasumiène brojeve èija je
# distribucija normalna. Argument x u rnorm funkciji je 
# broj nasumiènih brojeva koje æe funkcija vratiti. rnorm takoğer
# prihvaæa neobavezen argumente za srednju vrijednost (mean) i 
# standardnu devijaciju (sd):
rm(list=ls())
y <- rnorm(1000)
y
hist(y)

y1 <- rnorm(1000, mean=3)
y1
hist(y1)

y2 <- rnorm(1000, mean=3, sd=3)
y2
hist(y2)


#############
# 8. excerise
#############
rm(list=ls())
par(mfrow=c(2,1), mar=c(2,2,2,1))
#plot1
mean1 <- 0
sd1 <- 2
a <- mean1+sd1
b <- mean1-sd1
x <- seq(-4,4,length=200)*sd1+mean1
y <- dnorm(x,mean1,sd1)
plot(x, y, type="l", bty="l")
i <- x >= b & x <= a
polx <- c(b,x[i],a)
poly <- c(0,y[i],0)
polygon(polx,poly, col="red")

#plot 2
mean2 <- 3
sd2 <- 1
a2 <- mean2+sd2
b2 <- mean2-sd2
x2 <- seq(-4,4,length=200)*sd2+mean2
y2 <- dnorm(x2,mean2,sd2)
plot(x2, y2, type="l", bty="l")
i2 <- x2 >= b2 & x2 <= a2
polx2 <- c(b2,x2[i2],a2)
poly2 <- c(0,y2[i2],0)
polygon(polx2,poly2, col="green")

#############
# 9. excerise
#############
rm(list=ls())
x <- rnorm (100000, mean=1, sd=2)
i <- x > 5
y <- x[i]
l <- length(y)
l #broj brojeva veæih od 5
# Raèunanje koliko brojeva od 100000 brojeva je veæe od 5
# pomoæu funkcije pnorm s argumentom lower.tail=FALSE.

y1 <- pnorm (5, mean=1, sd=2, lower.tail=FALSE)
y1
# y1 oznaèava vjerojatnost da je broj dobiven normalnom
# distribucijom veæi od 5 pri mean=1 i sd=2 

y2 <- round (y1*100000)
y2
# mnoenjem vjerojatnosti y1 sa 100000 dobijemo broj brojeva 
# veæih od  5 

#plot1
par(mfrow = c(1,1))
mean1 <- 1
sd1 <- 2
x <- seq(-4,4,length=2000)*sd1+mean1
y <- dnorm(x,mean1,sd1)
plot(x, y, type="l", bty="l")

i <- x <= mean1-sd1
polx <- c(-4, x[i], mean1-sd1)
poly <- c(0, y[i], 0)
polygon(polx,poly, col="red")

j <- x >= mean1+sd1
polx2 <- c(mean1+sd1, x[j], 7)
poly2 <- c(0, y[j], 0)
polygon(polx2,poly2, col="red")

##############
# 10. excerise
##############
par(mfrow=c(2,1), mar=c(2,2,2,1))

#plot 1
#CDF
mean <- 0
sd <- 1
x <- seq(-4, 4, length.out = 1000) * sd + mean
y1 <- pnorm(x, mean, sd)
plot(x, y1, type = "l", bty = "l", col="red", ylim=c(0,1))
#f(x)=x
lines(x, x, col = "blue")
#qnorm 
y2 <- qnorm(x, mean, sd)
lines(x, y2, col="green")

rm(list=ls())

#plot 2
#CDF
mean <- 2
sd <- 3
x <- seq(-4, 4, length.out = 1000) * sd + mean
y1 <- pnorm(x, mean, sd)
plot(x, y1, type = "l", bty = "l", col="red", ylim=c(0,1))
#f(x)=x
lines(x, x, col = "blue")
#qnorm 
y2 <- qnorm(x, mean, sd)
lines(x, y2, col="green")

#qnorm i pnorm su meğusobno inverzne funkcije
x <- 0.5
a <- qnorm (x)
b <- pnorm (a)
all.equal(x,b)

qnorm(pnorm(3))
pnorm(qnorm(0.4))

#svejedno je napišemo li kompoziciju qnorm(pnorm(x)) ili 
#pnorm(qnorm(x)), u oba sluèaja vraæa se argument s kojim se
#ulazi u kompoziciju 

##############
# 11. excerise
##############
rm(list=ls()) 
par(mfrow = c(1,1))
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

##############
# 12. excerise
##############
rm(list=ls())

# uzorak 1
s1 <- rpois (100000, lambda = 6.5) 
# mean je jednak lambda

sd1 <- sd(s1)
sd1 
# standardna devijacija uzorka 1 (mean = 6.5)

var1 <- var(s1)
var1 
# varijanca uzorka 1 (mean = 6.5)

dif1 <- abs(var1-6.5)
dif1
# vrijednost varijance je vrlo bliska srednjoj vrijednosti uzorka 

# uzorak 2
s2 <- rpois (100000, lambda = 4.7) 
var2 <- var(s2)
dif2 <- abs(var2-4.7)
dif2

# uzorak 3
s3 <- rpois (100000, lambda = 3.2) 
var3 <- var(s3)
dif3 <- abs(var3-3.2)
dif3

# uzorak 4
s4 <- rpois (100000, lambda = 4.11) 
var4 <- var(s4)
dif4 <- abs(var4-4.11)
dif4

# uzorak 5
s5 <- rpois (100000, lambda = 1.23) 
var5 <- var(s5)
dif5 <- abs(var5-1.23)
dif5

# vrijednosti varijance svih uzoraka su vrlo bliske srednjoj 
# vrijednosti uzorka, razlika te dvije vrijednosti je vrlo mala

##############
# 13. excerise
##############
rm(list=ls())
set.seed (2)
a <- rnorm(5)
a

# set.seed je funkcija koja postavlja "seed" za generator 
# nasumiènih brojeva (RNG) u R-u. RNG u R-u je ustvari pseudo-RGN
# koji generira brojeve algoritmom koji za isti "seed" daje 
# uvijek istu sekvencu brojeva. Tako su nasumièno generirani 
# brojevi pod istim "seedom" jednaki i ostaju jednaki bez obzira 
# na duljinu sekvence. 

# Korištenje set.seed funkcije je korisno kada elimo nasimiènu 
# sekvencu brojeva koja æe biti ponovljiva svaki put.  

set.seed (2)
b <- rnorm(5)
b
a

# a i b sadre jednake brojeve jer su generirani s istim "seedom"

set.seed (2)
d <- rnorm(7)
d
a

# d sadri istu sekvencu prvih 5 brojeva kao a i b, a sekvenca se
# nastavlja i dalje jer je d dui od a i b

set.seed(1)
e <- rnorm(5)
e
a

# e sadri razlièite brojeve od a, b i d jer je generiran s
# razlièitim seedom 

##############
# 14. excerise
##############
# Faktor je varijabla u koju su vrijednosti spremljene kao 
# niz integera sa pridruenim setom character vrijednosti
# koje se koriste kad se faktor prikazuje (levels). 
# U faktoru je svaka jedinstvena vrijednost spremljena samo
# jednom (kategorija), a niz integera pokazuje kako su jedinstvene
# vrijednosti poredane u faktoru. 

rm(list=ls())
data = c("a", "b", "c", "d", "e", "f", "a", "b", "g", "e", "a")
f.data = factor(data)
f.data

data2 = c("a", "b", "c", "d", "g", "e", "a")
f.data2 = factor(data2)
f.data2

# U vektor su vrijednosti spremljene po redu kako ih unesemo, 
# nema kategorija nego je svaka vrijednost spremljena baš onako
# kako je unesemo u vektor. Vektor je niz podataka indeksiranih
# po poziciji u jednoj dimenziji. Osim toga, na vektorima moemo
# koristiti aritmetièke funkcije, a na faktorima ne (jer su
# faktori kategorizirani). 

##############
# 15. excerise
##############
# my.rnorm
rm(list=ls())
my.rnorm <- function (n, mean, sd){
  set.seed(2)
  b <- runif(n)
  d <- qnorm(b, mean, sd)
  return(d)
}
my.rnorm(3, 0, 1)

# test
set.seed(2)
rnorm(3, 0, 1)
# my.rnorm vraæa neke iste brojeve kao i rnorm za isti n 
# (uz isti seed), ali neki brojevi nisu isti, ne znam zašto :)


# my.rpois
rm(list=ls())
my.rpois <- function (n, l){
  set.seed(3)
  b <- runif(n)
  d <- qpois(b, l)
  return(d)
}
my.rpois(13, 3)

# test
set.seed(3)
rpois(13, 3)
# my.rpois vraæa iste brojeve kao i rpois za isti n i lamda 
# (uz isti seed)

# my.runif
rm(list=ls())
my.runif <- function (n, min, max){
  a <- seq(min, max, length = 20000000)
  set.seed(3)
  b <- sample(a, n)
  return(b)
}
my.runif(3, 1, 5)

# test
set.seed(3)
runif(3, 1, 5)
# my.runif vraæa iste brojeve kao i runif za isti n, min i max 
# (uz isti seed)

##############
# 16. excerise
##############
rm(list=ls())
f.letters <- factor(letters)
f.letters
levels(f.letters) = c("0", "1", "2", "d", "e", "f", "g", "h", "i", 
                      "j", "k", "l", "m", "n", "o", "p", "q", "r", 
                      "s", "t", "u", "v", "w", "x", "y", "z")
f.letters
str(f.letters)

# str prikazuje strukturu varijable f.letters. Pokazuje broj 
# levela u factoru (26), pokazuje koji su to leveli kao niz 
# charactera i takoğer pokazuje niz brojeva. Kako bi bolje vidio
# što taj niz brojeva znaèi napravio sam donji primjer s manje 
# slova:

letters2 <- c("a", "b", "c", "d", "a", "a", "b", "d")
f.letters2 <- factor(letters2)
f.letters2
str(f.letters2)

# Iz ovog primjera zakljuèujem da leveli u faktoru predstavljaju
# kategorije vrijednosti s kojima ulazimo u faktor. U vektoru 
# letters2 imamo 3 puta slovo "a", a kada vektor pretvorimo u faktor
# pod levels se a pojavljuje samo jednom. Niz brojeva nakon popisa
# levela pokazuje kako su pojedine vrijednosti poredane u faktoru: 
# svaki broj predstavlja broj mjesta u levels na kojem se nalazi
# vrijednost za to mjesto. Npr. ako je levels: "a", "b", "c",
# onda struktura 1 1 2 1 3 znaèi da je faktor: a a b a c:

letters3 <- c("a", "a", "b", "a", "c")
f.letters3 <- factor(letters3)
str(f.letters3)
f.letters3

##############
# 17. excerise
##############
rm(list=ls())
class(iris)
iris

# class objekta iris je data.frame


a <- nrow(iris)
a
b <- ncol(iris)
b

# broj redova data.frame-a se moe saznati funkcijom nrow, a broj 
# stupaca faunkcijom ncol. Broj redova je 150, a broj stupaca 5. 


sapply(iris, class)

# class varijabli unutar data.frame-a moemo saznati tako da
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

# virginica je data.set koji sadri podatke iz iris samo za 
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

##############
# 18. excerise
##############
rm(list=ls())
m1 <- matrix(1:30, nrow = 5, ncol = 6, byrow = TRUE)
m1

m2 <- matrix(1:6, nrow = 2, ncol = 3)
m2

m1 [4:5, 4:6] <- m2
m1

mx1 <- m1
mx1 <- mx1[, -2]
mx1

##############
# 19. excerise
##############
# Najveæa razlika izmeğu data frame i matrice je što data frame moe
# sadravati razlièite vrste podataka, npr. u jednom stupcu mogu biti
# brojevi, u drugom characteri, u treæem faktori itd. 
# Matrice mogu sadravati samo jedan tip podataka, tj. svi redovi i stupci
# matrice moraju biti imati isti class.

##############
# 20. excerise
##############
rm(list=ls())
v1 <- seq(10, 60, length = 6)
v1

m1 <- matrix(1:30, nrow = 5, ncol = 6, byrow = TRUE)
m1

umn1 <- m1 * v1
umn1

# m1 * v1 - mnoe se elementi matrice i vektora s obzirom na njihovu poziciju
# Mnoenje ide po redovima: element prvog reda u prvom stupcu matrice se 
# mnoi s prvim elementom vektora, element drugog reda u prvom stupcu s 
# drugim elementom vektora itd. S obzirom da matrica ima 5 redova, a vektor 
# 6 elemenata, šesti element vektora se mnoi s prvim redom drugog stupca.
# Zatim se opet prvi element vektora mnoi s elementom drugog reda i drugog 
# stupca matrice itd. sve dok se svi elementi matrice ne pomnoe s vektorom
# (recikliranje)

umn2 <- m1 %*% v1
umn2
m1[2,]*v1
sum(m1[2,]*v1)


# m1 %*% v1 - promatra vektor kao matricu s jednim stupcem i 6 redova,
# zatim mnoi matricu m1 s matricom v1:

#      [,1] [,2] [,3] [,4] [,5] [,6]        [,1]
# [1,]    1    2    3    4    5    6   [1,]   10
# [2,]    7    8    9   10   11   12   [2,]   20
# [3,]   13   14   15   16   17   18 x [3,]   30
# [4,]   19   20   21   22   23   24   [4,]   40
# [5,]   25   26   27   28   29   30   [5,]   50
#                                      [6,]   60
# 
# Mnoe se redovi sa stupcem na ovaj naèin:
# 1*10 + 2*20 + 3*30 + 4*40 + 5*50 + 6*60 = 910
# 7*10 + 8*20 + 9*30 + 10*40 + 11*50 + 12*60 = 2170 
# tj. prvi element stupca s prvim elementom reda, drugi element s drugim itd.
# rezultat je matrica s jednim stupcem i 5 redova. 

umn3 <- v1 %*% m1

# nije moguæe obaviti operaciju. 
# v1 se promatra kao matrica s jednim redom i 6 stupaca (6 elemenata), 
# a svaki stupac matrice m1 ima 5 redova (5 elemenata). Algoritam treba 
# pomnoiti red matrice v1 sa svakim stupcem matrice m1, ali su razlièite
# duljine (non-conformable) pa to nije moguæe. Kad bi npr. vektor v1 imao 
# 5 elemenata tada bi v1 %*% m1 bilo moguæe:

v1 <- seq(10, 50, length=5)
v1
umn3 <- v1 %*% m1
umn3

##############
# 21. excerise
##############
rm(list=ls())
v1 <- seq(10, 60, length = 6)
v1

m1 <- matrix(1:30, nrow = 5, ncol = 6, byrow = TRUE)
m1

zbr1 <- m1 + v1
zbr1

# m1 * v1 - zbrajaju  se elementi matrice i vektora s obzirom na njihovu poziciju
# Zbrajanje ide po redovima: element prvog reda u prvom stupcu matrice se 
# zbraja s prvim elementom vektora, element drugog reda u prvom stupcu s 
# drugim elementom vektora itd. S obzirom da matrica ima 5 redova, a vektor 
# 6 elemenata, šesti element vektora se zbraja s prvim redom drugog stupca.
# Zatim se opet prvi element vektora zbraja s elementom drugog reda i drugog 
# stupca matrice itd. sve dok se svi elementi matrice ne zbroje s vektorom
# (recikliranje)

##############
# 22. excerise
##############
rm(list=ls())
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
# mnoenjem mt1 %*% mt2 m i p, tada su dimenzije matrice
# mt2 n i p. 

n <- 4
m <- 5
p <- 7
mt1 <- matrix(1:100, nrow = m, ncol = n)
mt2 <- matrix(1:100, nrow = n, ncol = p )
mt3 <- mt1 %*% mt2
a <- dim(mt3)
a[1] == m & a[2] == p


# Generalna pravila kod mnoenja:
# - m1 %*% m2 = m3: 
#   - dimenzije m1 su m i n
#   - tada su dimenzije m2 n i p
#   - tada su dimenzije m3 m i p

# - m2 %*% m1 = m3
#   - dimenzije m2 n i p
#   - tada su m1 su p i m
#   - tada su dimenzije m3 n i p

##############
# 23. excerise
##############
rm(list=ls())
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

##############
# 24. excerise
##############
# Lista je vektor :)
rm(list=ls())
lst <- list(vc1 = c(1, 3, 5), vc2 = c(2, 4, 6), vc3 = c(10, 20, 30, 40, 50))
is.vector(lst)

# U listu se moe spremiti više tipova podataka (pa èak i druge liste), 
# dok se u pojedinaène vektore moe spremiti samo jedan tip podataka 
# (numeric, character, logical)

##############
# 25. excerise
##############
rm(list=ls())
iris
a <- iris[["Sepal.length"]]
a
b <- iris[1]
d <- iris[[1]]
d

class(a)
class(b)
class(d)

# iris[["Sepal.length"]] vraæa NULL
# iris[1] vraæa data.frame s vrijednostima iz prvog stupca data frame iris
# iris[[1]] vraæa numerièke vrijednosti iz prvog stupca iris 

is.list(iris)
# data.frame je vrsta liste

##############
# 26. excerise
##############
# diag zamjenjuje dijagonalu matrice, a moe se koristiti i za dobivanje
# vrijednosti dijagonale matrice (diag(m), m je matrix)

rm(list=ls())
mm1 <- matrix (1:9, nrow = 3, ncol = 3, byrow = TRUE)
mm1

diag (mm1) <- c(10, 20, 30)
mm1

# syntax je èudan, kao da je diag varijabla unutar matrice kojoj
# pridruujemo vrijednost.
# Google kae: 
# diag(object) <- value je skraæeno od object <- "diag<-"(object, value)
# tj. diag (mm1) <- c(10, 20, 30) je"mm1 <- "diag<-"(mm1, c(10, 20, 30)). 

rm(list=ls())
mm1 <- matrix (1:9, nrow = 3, ncol = 3, byrow = TRUE)
mm1 <- "diag<-"(mm1, c(10, 20, 30))
mm1

# Sve replacements funkcije u R-u su napravljene tako da mi pozivamo skraæeni
# oblik, a R ustvari preimenuje naš skraæeni oblik u dui oblik i zove 
# odreğenu replacements funkciju. Opæenito:
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
