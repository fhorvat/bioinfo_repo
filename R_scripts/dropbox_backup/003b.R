# CpG otoci su regije u DNA s visokim sadržajem 
# C+G nukleotida u kojima se èesto pojavljuju CpG
# dinukleotidi (s obzirom na ukupnu kolièinu DNA). 
# Inaèe je DNA kralješnjaka siromašna dinukleotidom
# CpG. 
# Definicija iz rada: CpG otok je DNA sekvenca u kojoj je:
# - pomièni prosjek %G+C baza veæi od 50
# - pomièni prosjek omjera primjeæenih CpG/oèekivanih CpG
# (Obs/Exp CpG) je veæi od 0.6

# Obs/Exp CpG = ((number of CpG)/(number of C * number of G)) * N
# gdje je N ukupan broj nukleotida u analiziranoj sekvenci

# U radu promatraju sekvencu pomoæu "prozora" od 100 bp (N=100)
# koji se pomièe jednu po jednu bazu. 

CpG <- scan("141107_CpG", what = "character", sep = "\n")
a <- grep(">", CpG)
CpG <- CpG[-a]
CpG <- paste(CpG, sep = "", collapse = "")


# Izraèun pomiènog zbroja CpG dinukleotida. Zbrajam CpG dinukleodite
# u bazama 1:100, pa zatim 2:101, pa 3:102 itd. 

# Prvo pomoæu vlastite funkcije only.CG sve CG dinukleotide
# zamjenjujem s 12, a onda sve 2 s 0. To radim da bi saèuvao
# poredak nukeotida (zato stavljam dvoznamenkasti broj), a zatim
# zamjenjujem 2 s 0 da bi mogao "prebrojati" CG zbrajanjem jedinica.
# Nakon funkcije je praktièki svaki CpG dinukleotid zamijenjen
# sa 1, ali tako da iza slijedi 0 da bi poredak nukleotida ostao
# saèuvan. Ostali nukleotidi su zamijenjeni s 0. 

# Zatim cijepam character vektor duljine 1 na vektor u kojem
# je svaki character odvojen. Na kraju na dobivenim brojevima
# funkcijom move.sum pomièno zbrajam CpG dinukleotide. 

# Funkcija pomiènog zbroja move.sum:
# pomièe se po argumentu x za 1 uzimajuæi 100 èlanova i 
# zbraja ih. Koristi funkciju filter koja ustvari radi na 
# time-series  podacima, ali ju koristim jer nisam našao neku 
# drugu funkciju koja je "kotrljajuæa", a sve petlje koje sam
# probao su prespore na ovako puno podataka.
# Funkcija vraæa objekt klase Time-Series pa sve pretvaram u 
# brojeve.

# Cijelo vrijeme koristim istu varijablu (maštovitog naziva xyz)
# kako mi se ne bi nakupljale velike varijable i zauzimale memoriju.
# Ionako sve odmah stavljam u data.frame. 

only.CG <- function(x){
  y <- gsub("CG", 12, x)
  y <- gsub(2, 0, y)
  y <- gsub("[GATC]", 0, y)
  return(y)
}

move.sum <- function(x, n = 100){
  filter(x, rep(1, n), sides = 1)
}

xyz <- only.CG(CpG)
xyz <- unlist(strsplit(xyz, split = ""))
xyz <- as.numeric(move.sum(xyz))
CpG.df <- data.frame(xyz)
colnames(CpG.df) <- "Number of CpG"

# Zbrajanje C i G nukleotida radim malo drukèije, 
# pretvaranje nukleotida u brojeve radim pomoæu faktora. 
# Mislim da je ovo brži naèin nego da radim analogno ovom iznad
# (gsub).

# Prije svega vektor duljine 1 sa svim bazama pretvaram
# u vektor duljine jednakoj duljini sekvence u kojoj je svaka
# baza zasebni èlan vektora.

CpG <- unlist(strsplit(CpG, split = ""))

xyz <- factor(CpG) # levels(CG.mov.avg) = A C G T
levels(xyz) <- c(0, 1, 0, 0) 
xyz <- as.character(xyz)
xyz <- as.numeric(move.sum(xyz))
CpG.df <- cbind(CpG.df, xyz)
colnames(CpG.df) <- c("Number of CpG", "Number of C")

xyz <- factor(CpG) # levels(CG.mov.avg) = A C G T
levels(xyz) <- c(0, 0, 1, 0) 
xyz <- as.character(xyz)
xyz <- as.numeric(move.sum(xyz))
CpG.df <- cbind(CpG.df, xyz)
colnames(CpG.df) <- c("Number of CpG", "Number of C", "Number of G")

# U kolonu Obs/Exp stavljam rezultate izraèuna pomiènog prosjeka
# primjeæenih CpG/oèekivanih CpG (Obs/Exp CpG):
CpG.df$"Obs/Exp" <- (CpG.df$"Number of CpG"/
                       (CpG.df$"Number of C"
                        *CpG.df$"Number of G"))*100

# Izraèun pomiènog prosjeka C+G. Radi slièno kao i izraèun 
# pomiènog zbrajanja samo što funkcija raèuna pomièni prosjek;
# uzima prvih èlanove od 1:100 i raèuna njihov prosjek, zatim
# 2:101, pa 3:102 itd. Isto kao i move.sum koristi funkciju 
# factor (i iz istog razloga). Levels faktora su pomnoženi sa
# 100 kako bi dobio postotak. 

move.average <- function(x, n = 100){
  filter(x, rep(1/n, n), sides = 1) 
}

xyz <- factor(CpG) # levels(CG.mov.avg) = A C G T
levels(xyz) <- c(0, 100, 100, 0) 
xyz <- as.character(xyz)
xyz <- as.numeric(move.average(xyz))
CpG.df <- cbind(CpG.df, xyz)
colnames(CpG.df) <- c("Number of CpG", "Number of C", 
                      "Number of G", "Obs/Exp", "%C+G")
