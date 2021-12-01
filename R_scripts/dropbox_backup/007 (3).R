# Funkcije Xnorm koriste se za operacije vezane uz
# normalnu distribuciju. 

# Normalna distribucija je kontinuirana distribucija 
# vjerojatnosti koja govori kolika je vjerojatnost da 
# neka primjeæena realna vrijednost padne izmeðu bilo 
# koje dvije realne granice. 

# dnorm je funkcija koja uzima argument x i vraæa visinu
# distribucije vjerojatnosti za svaku toèku unutar x. 
# Ako se koristi samo kao dnorm(x) uzima srednju vrijednost 0 
# i standardnu devijaciju 1:
x <- seq(-4,4,length=200)
y <- dnorm (x)
plot(x,y,type="l",col="red")

# Moguæe je i definirati proizvoljnu srednju vrijednost i standardnu
# devijaciju: dnorm (x, mean=, sd=). Primjer:
x <- seq(-10,10,length=200)
y <- dnorm (x, mean=2, sd=2)
plot(x,y,type="l",col="red", xlim=c(-10,10))

# Ukoliko standardna devijacija (sd) ostaje nepromjenjena, 
# poveæavanjem srednje vrijednosti (mean) dolazi do horizontalnog 
# pomicanja krivulje udesno, a smanjenje mean pomièe krivulju 
# horizontalno ulijevo:
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
x <- seq(-15,15,length=200)
y1 <- pnorm (x, mean=2, sd=4)
y2 <- pnorm (x)
plot (x,y1,type="l",col="red", xlim=c(-15,15), ylim=c(0,1))
lines(x,y2,col="green")

# Ako želimo da funkcija daje vjerojatnost da je nasumièan broj 
# dobiven normalnom distribucijom VEÆI od x možemo korisiti 
# argument lower.tail=FALSE:
x <- seq(-15,15,length=200)
y1 <- pnorm (x, mean=2, sd=4, lower.tail=FALSE)
y2 <- pnorm (x, lower.tail=FALSE)
plot (x,y1,type="l",col="red", xlim=c(-15,15), ylim=c(0,1))
lines(x,y2,col="green")

# qnorm je funkcija inverzna funkciji pnorm: vrijednost argumenta
# x je vjerojatnost, a funkcija vraæa broj èija kumulativna 
# distribucija odgovara toj vjerojatnosti x. qnorm takoðer prihvaæa
# neobavezne argumente mean i sd:
x <- seq(0,1,length.out=100)
y <- qnorm(x)
y2 <- qnorm(x, mean=2, sd=3)
plot(x,y, type="l",col="red", ylim=c(-10,10))
lines (x,y2, col='yellow')

# rnorm je funkcija koja generira nasumiène brojeve èija je
# distribucija normalna. Argument x u rnorm funkciji je 
# broj nasumiènih brojeva koje æe funkcija vratiti. rnorm takoðer
# prihvaæa neobavezen argumente za srednju vrijednost (mean) i 
# standardnu devijaciju (sd):

y <- rnorm(x)
y
hist(y)

y1 <- rnorm(x,mean=3)
y1
hist(y1)

y2 <- rnorm(x,mean=3,sd=3)
y2
hist(y2)




