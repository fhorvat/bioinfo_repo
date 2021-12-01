#############
# 1. excerise
#############
# Normal distribution:
par(mfrow = c(1, 1))
rm(list=ls())

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
    rnd.x[i] <- list(matrix(rnorm(n[i]*10000), nrow = n[i], ncol = 10000, 
                            byrow = FALSE))
    mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
  }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

# plot distribution of means:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))

# standard deviation
# Standardna devijacija je prikazana u data.frame-u na zadnjem mjestu 
# u listi koju vraæa funkcija Nor.dist. Što je veæi n, standardna
# devijacija je manja (obrnuto proporcionalno)


# Uniformna distribucija 
par(mfrow = c(1, 1))
rm(list=ls())

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
    rnd.x[i] <- list(matrix(runif(n[i]*10000), nrow = n[i], ncol = 10000, 
                            byrow = FALSE))
    mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
  }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

# plot distribution:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))


# Poissonova distribucija 
par(mfrow = c(1, 1))
rm(list=ls())

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
    rnd.x[i] <- list(matrix(rpois(n[i]*10000, l = 2), nrow = n[i], 
                            ncol = 10000, byrow = FALSE))
    mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
  }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

# plot distribution:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))


# Moja distribucija 
par(mfrow = c(1, 1))
rm(list=ls())

# vraæa n brojeva po distribuciji, brojevi idu [-2*pi, 2*pi], 
# a vjerojatnosti od [(sin(-2*pi))^2, (sin(2*pi))^2]
My.random.gen <- function(n){
  rnd.gen <- seq(-2*pi, 2*pi, length.out = 10000)
  rnd.gen.probability <- (sin(rnd.gen))^2
  xyz <- sample(rnd.gen, n, prob = rnd.gen.probability, replace = TRUE)
  return(xyz)
}

# plot moje distribucije, napravio sam ga u posebnoj funkciji da se ne
# bi kasnije kad pozovem funkciju unutar funkcije Distribution i puta
# iscrtavao plot:
My.random.gen.plot <- function(n){
  rnd.gen <- seq(-2*pi, 2*pi, length.out = 10000)
  rnd.gen.probability <- (sin(rnd.gen))^2
  plot(rnd.gen, rnd.gen.probability, type = "l")
}
My.random.gen.plot(1)

Distribution <- function(n){
  rnd.x <- list()
  mean.rnd.x <- list()
  for (i in 1:length(n)){
    rnd.x[i] <- list(matrix(My.random.gen(n[i]*10000), nrow = n[i], 
                            ncol = 10000, byrow = FALSE))
    mean.rnd.x[i] <- list(colMeans(rnd.x[[i]]))
  }
  stn.dev <- 0
  for (j in 1:length(n)){
    stn.dev[j] <- sd(unlist(mean.rnd.x[j]))
  }
  mean.rnd.x$standard.deviation <- data.frame(n, stn.dev)
  plot(stn.dev, n, type = 'l')
  return(mean.rnd.x)
}
set.seed(1)
Distribution(c(5, 50, 500, 5000))

# plot distribution:
a <- Distribution(c(5, 50, 500, 5000))
par(mfrow = c(2, 2))
hist(a[[1]])
hist(a[[2]])
hist(a[[3]])
hist(a[[4]])
par(mfrow = c(1, 1))

#############
# 2. excerise
#############
rm(list=ls())
# a)
chilled.Mc1 <- CO2[CO2$Plant == "Mc1" & CO2$Treatment == "chilled", ]
chilled.Mc1

# b)
chilled.Q.less250CO2 <- CO2[CO2$Type == "Quebec" & CO2$Treatment == "chilled"
                            & CO2$conc < 250, ]
chilled.Q.less250CO2

# c)
chilled.Q.less250CO2.odd <- CO2[c(TRUE, FALSE) & CO2$Type == "Quebec" 
                                & CO2$Treatment == "chilled"
                                & CO2$conc < 250,  ]
chilled.Q.less250CO2.odd

# d) 
more350CO2.more35up <- CO2[CO2$conc > 350 & CO2$uptake > 35, ] 
more350CO2.more35up

# e) 
more350CO2.more35up.plant.type <- CO2[CO2$conc > 350 & CO2$uptake > 35,
                                      c(1, 2)]
more350CO2.more35up.plant.type

# f)
# a)
sub.chilled.Mc1 <- subset(CO2, Plant == "Mc1" & Treatment == "chilled")
sub.chilled.Mc1

# b)
sub.chilled.Q.less250CO2 <- subset (CO2, Type == "Quebec" 
                                    & Treatment == "chilled" & conc < 250)
sub.chilled.Q.less250CO2 

# c) 
sub.chilled.Q.less250CO2.odd <- subset(CO2, c(TRUE, FALSE) & Type == "Quebec" 
                                       & Treatment == "chilled" & conc < 250)
sub.chilled.Q.less250CO2.odd 

# d) 
sub.more350CO2.more35up <- subset(CO2, conc > 350 & uptake > 35)
sub.more350CO2.more35up

# e)
sub.more350CO2.more35up.plant.type <- subset(CO2, conc > 350 & uptake > 35, 
                                             c(Plant, Type))
sub.more350CO2.more35up.plant.type

#############
# 3. excerise
#############
rm(list=ls())
x <- sum(rep(c(T, T, F), times = 5))
x

a <- rep(c(T, T, F), times = 5)
a

b <- sum(c(T, T, F))
b

# - rezultat sum(rep(c(T, T, F), times = 5)) je 10

# OBAŠNJENJE:

# - rep sekvencu T T F replicira 5 puta

# - nakon rep(c(T, T, F), times = 5) imamo T T F sekvencu 5 puta tj. 
# 10 puta T i 5 puta F

# - sum radi na numeric argumentima

# - kad pomoæu sum pokušamo zbrojiti logic argumente R ih pretvori u 
# numeric tako da je TRUE = 1, a FALSE = 0

# - kako je T = 1, a F = 0, kada sekvencu (T T F) x 5 koja ima 10xT i 5xF
# zbrojimo pomoæu sum imamo 10*1 + 5*0 što daje 10


#############
# 4. excerise
#############
rm(list=ls())
my.data.read.table <- read.table("gene_with_protein_product.txt", header = T,
                                 sep =  "\t", fill = TRUE, quote = "")

my.data.read.csv <- read.csv("gene_with_protein_product.txt", header = T, 
                             sep =  "\t")

my.data.read.delim2 <- read.delim2("gene_with_protein_product.txt")

# read.csv i read.delim kao oznaku za decimalni znak èitaju
# ".", a read.csv2 i read.delim2 ","

# stringsAsFactors je logical argument koji kaže R-u treba li
# stringove iz filea kojeg uèitavamo pretvoriti u faktore
# (ako je TRUE) ili u charactere (ako je FALSE).
# Pretpostavljena vrijednost je TRUE. 

my.data.Fact <- read.table("gene_with_protein_product.txt", header = T,
                           sep =  "\t", fill = TRUE, quote = "")
a <- sapply(my.data.Fact, class)
a

my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)
b <- sapply(my.data.Char, class)
b

#############
# 5. excerise
#############
a <- sum(nchar(my.data.Char$Approved.Symbol)<5)
a 

#############
# 6. excerise
#############
a <- grep("^HGNC", my.data.Char$HGNC.ID, value = T) == my.data.Char$HGNC.ID
a

# a je logièka varijabla u kojoj su T i F vrijednosti usporedbe HGNC.ID
# stupca my.data.Char data.frame-a i character varijable dobivene grep 
# funkcijom koja iz tog istog stupca uzima sve vrijednosti koje poèinju s
# HGNC

all(a)
# all vraæa TRUE ako su svi èlanovi TRUE, tj. vraæa FALSE ako je 
# barem jedan èlan FALSE

# all(a) æe vratiti TRUE ako sve vrijednosti stupca $HGNC.ID poèinju
# s HGNC

#############
# 7. excerise
#############
a <- grep("\\d$", my.data.Char$Approved.Symbol, value = T)
a
# a se sastoji do svih approved symbols koji na zadnje mjestu imaju
# znamenku


b <- grep("\\d.$", my.data.Char$Approved.Symbol, value = T)
b
# b se sastoji do svih approved symbols koji na predzadnjem mjestu 
# imaju znamenku

# value je parametar funkcije grep koji ako je TRUE vraæa vrijednost
# koja odgovara zadanom patternu. Ako je FALSE vraæa
# mjesto (lokaciju, index) na kojoj se nalazi vrijednost koja 
# odgovara patternu. Pretpostavljena vrijednost je FALSE. Primjer:

xyz <- c("plavo", "zeleno", "crveno", "žuto", "bijelo")
x <- grep("l", xyz, value = TRUE)
y <- grep("l", xyz, value = FALSE)
x
y

# x sadrži sve rijeèi iz character vektora xyz koje sadrže slovo l
# y sadrži broj mjesta (smještaj, index) u character vektoru xyz 
# na kojem se nalazi rijeè koja sadrži slovo l 

#############
# 8. excerise
#############
xyz <- c("plavo", "zeleno", "crveno", "žuto", "bijelo")

# Kao što sam objasnio u prošlom zadatku, grep vraæa vrijednost
# iz character vektora koja odgovara zadnom patternu ili vraæa 
# smještaj te vrijednosti u vektoru (ovisno o argumentu value)

# grepl u character vektoru provjerava jednu po jednu vrijednost
# za zadani pattern i vraæa TRUE ili FALSE za svaku vrijednost 
# vektora s obzirom na to sadrži li zadani pattern
# (TRUE ako sadrži, FALSE ako ne sadrži). Primjer:

x <- grep("v", xyz, value = TRUE)
y <- grepl("v", xyz)
x
y

# U primjeru provjeravamo koja rijeè iz vektora xyz sadrži slovo "v". 
# x je dobiven grep funkcijom i sadrži sve rijeèi koje u sebi
# imaju slovo "v". 
# y je dobiven grepl funkcijom i sadrži logièke vrijednosti po redu
# za svaku rijeè iz vektora xyz s obzirom na to sadrži li slovo "v" ili ne.
# Npr. na prvom mjestu u xyz je vrijednost "plavo" koja sadrži "v" 
# pa je prva vrijednost u logièkom vektoru y TRUE. Na drugom mjestu
# u xyz je "zeleno" što ne sadrži "l" pa je zato na drugom mjestu u 
# logièkom vektoru y FALSE itd. za sve vrijednosti iz vektora xyz.

#############
# 9. excerise
#############
my.data.Char$my.hgnc <- substr(my.data.Char$HGNC.ID, 6, 
                               nchar(my.data.Char$HGNC.ID))

##############
# 10. excerise
##############
a <- my.data.Char[c(grep("RNASELI", my.data.Char$Approved.Symbol),
                    grep("RNASELI", my.data.Char$Previous.Symbols)), ]
a

# a sadrži redove iz data.frame-a u èijim stupcima Approved.Symbol
# i Previous.Symbols se pojavljuje gen RNASELI

##############
# 11. excerise
##############
my.data.Char$my.approved.symbol <- sub("CYP", "XXX", 
                                       my.data.Char$Approved.Symbol)

# provjera:
x <- my.data.Char[grep("XXX", my.data.Char$my.approved.symbol), ]
x

# sub zamjenjuje pattern u stringu zadanom zamjenom (replacement)
# prvi put kad se pattern pojavi u stringu

# gsub zamjenjuje pattern u stringu zadanom zamjenom svaki puta
# kad se pattern pojavi u stringu. Primjer:

xyz <- c("Ivan", "Ivana", "Anastazija", "Valentina", "Tea")
y <- sub("a", "@", xyz)
z <- gsub("a", "@", xyz)
y
z

# U y je prvo slovo "a" u imenima zamijeneno s @
# U z je svako slovo "a" u imenima zamijenjeno s @

##############
# 12. excerise
##############
x <- ifelse(my.data.Char$Previous.Symbols == "", 
            my.data.Char$Approved.Symbol, 
            paste(my.data.Char$Approved.Symbol, ",", 
                  my.data.Char$Previous.Symbols))
x <- sub(" ", "", x)
my.data.Char$all.symbols <- x

##############
# 13. excerise
##############
# a
"^\\d[A-z]{16}\\d+\\?"

#b 
"^[[:alnum:]]{3}[b-yB-Y]+\\?"

# c
"^([0-3]{3})\\1+[A-z]*"

# d
"\\d$"

# e
"^b.*\\d$"

##############
# 14. excerise
##############
write.table(my.data.Char, file = "my.data.Char1.txt",
            sep = ";")

write.table(my.data.Char, file = "my.data.Char2.txt",
            sep = "\t")

##############
# 15. excerise
##############
rm(list=ls())
FASTA <- function(location){
  fasta.data <- readLines (location, warn = FALSE)
  xyz <- grep("^>", fasta.data) # pozicije u fasta.data na kojima se nalaze imena gena 
  seq.a <- xyz + 1 # pozicije na kojima poèinju sekvence gena 
  seq.z <- xyz - 1 # pozicije na kojima završavaju sekvence gena
  seq.z <- seq.z[-1] # prva vrijednost je 0 pa ju odstranjujemo
  seq.z <- c(seq.z, length(fasta.data)) # dodamo na kraj poziciju na kojoj završava zadnja sekvenca
  
  sequence <- list()
  for (i in 1:length(seq.a)){
    sequence[i] <- paste(fasta.data[seq.a[i]:seq.z[i]], collapse="")
  }
  
  seq.name <- grep("^>", fasta.data, value = TRUE)
  seq.name <- gsub("^>", "", seq.name)
  sequence <- as.character(sequence)
  
  pr1 <- grepl("[^GATC]", sequence) 
  if (any(pr) == TRUE){
    warning("sekvenca nije u FASTA formatu")
  }
    
  dt.fasta <- data.frame(seq.name, sequence)
  return(dt.fasta)
}

fasta.data.frame <- FASTA("fasta.txt")
fasta.data.frame

# Moram priznati da mi je ideju kako napraviti ovaj zadatak dala 
# Dora (ovaj prvi dio u funkciji, kako dobiti pozicije na kojima 
# poèinju i završavaju DNA sekvence, varijable seq.a i seq.b). 
# Takoðer nisam stigao napraviti provjeru ima li file koji uèitavamo
# pravilni format (jer je sada nedjelja u 23:42) :)

##############
# 16. excerise
##############
rm(list=ls())
lst <- list(vct1 = 1:5, vct2 = 5:73, vct3 = seq(from = - pi, to = exp(1), length.out = 24))

mean.lapply <- lapply(lst, mean)
mean.lapply

mean.sapply <- sapply(lst, mean)
mean.sapply

# lapply vraæa rezultat u obliku liste na kojoj su vrijednosti
# srednjih vrijednosti za pojedinaène vektore s liste lst

# sapply je pojednostavljena funkcija koja vraæa vrijednosti
# srednjih vrijednosti za pojedinaène vektore s liste lst
# u obliku vektora s imenima vektora s liste lst

mean.lst.for <- list()
for (i in 1:length(lst)){
  mean.lst.for[i] <- mean(lst[[i]])
}
mean.lst.for

x <- 1
mean.lst.while <- list()
while (x < length(lst)+1){
  mean.lst.while[x] <- mean(lst[[x]])
  x <- x + 1
}
mean.lst.while


##############
# 17. excerise
##############
rm(list=ls())
mean.sepal.length <- by(iris$Sepal.Length, iris$Species, mean)
mean.sepal.length

##############
# 18. excerise
##############
rm(list=ls())
rivers
sum.for <- 0
for (i in 1:length(rivers)){
  if (rivers[i] > 650) {
    sum.for <- sum.for + rivers[i]
  }
  else{
    sum.for
  }
}
sum.for

sum.subset <- sum(subset(rivers, rivers > 650))
sum.subset

##############
# 19. excerise
##############
Monty.hall <- function(){
  nagrade <- c("koza", "koza", "auto")
  vrata <- sample(nagrade, 3)
  vrata <- vrata[-1]  # jer natjecatelj uvijek bira vrata 1
  x <- match("koza", vrata)  # iza kojih vrata se nalazi koza 
  vrata <- vrata[-x]  # voditelj otvara vrata iza kojih se nalazi koza
  # ostaju samo jedna vrata
  # to su vrata koja natjecatelj može izabrati kada mu je ponuðena zamjena
  return(vrata) 
}
# 100000 puta pokrenemo funkciju Monty.Hall, 
# rezultati se spremaju:

rezultat <- replicate(100000, Monty.hall()) 

# raèunam kolika je šansa da je natjecatelj pobjedio za oba sluèaja:

# natjecatelj se odluèuje na zamjenu i pobjeðuje 
# (što znaèi da je iza "zadnjih preostalih" vrata iz 
# funkcije bio auto):

zamjena.pobjeda <- rezultat[rezultat == "auto"]
length(zamjena.pobjeda) # broj pobjeda ako se natjecatelj odluèi na zamjenu
vjerojatnost.zamjena.pobjeda <- length(zamjena.pobjeda)/100000
vjerojatnost.zamjena.pobjeda


# natjecatelj ostaje pri prvom izboru i pobjeðuje
# (iza prvih vrata je auto, iza "zadnjih" vrata je koza):
prvi.izbor.pobjeda <- rezultat[rezultat == "koza"]
length(prvi.izbor.pobjeda) # broj pobjeda ako natjecatelj ne mijenja vrata
vjerojatnost.prvi.izbor.pobjeda <- length(prvi.izbor.pobjeda)/100000
vjerojatnost.prvi.izbor.pobjeda
