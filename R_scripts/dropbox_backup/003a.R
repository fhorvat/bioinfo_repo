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
# koji se pomièe jednu po jednu bazu pa i ja koristim 100 bp.

# Èitanje podataka u character vektor duljine 1 (jedan string):
CpG <- scan("141105_CpG_small", what = "character", sep = "\n")
CpG <- CpG[-1]
CpG <- paste(CpG, sep = "", collapse = "")

# Funkcija pomiènog zbroja move.sum:
# pomièe se po argumentu x za 1 uzimajuæi 100 èlanova i na 
# zbraja ih.
# Koristi funkciju filter koja ustvari radi na time-series 
# podacima, ali ju koristim jer nisam našao neku drugu funkciju 
# koja je "kotrljajuæa", a sve petlje koje sam probao su 
# prespore na ovoliko puno podataka. 
move.sum <- function(x, n = 100){
  filter(x, rep(1, n), sides = 1)
}

# Slijede funkcije koje pomièno zbrajaju CpG dinukleotide te 
# C i G nukleotid koji mi trebaju za izraèun pomiènog prosjeka
# omjera primjeæenih CpG/oèekivanih CpG (Obs/Exp CpG). Koristim
# gore definiranu funkciju move.sum koja radi na brojevima.
# Zato pomoæu gsub u funkcijama zamjenjuje odgovarajuæe nukleotide
# s 1 i 0 i onda ih mogu pomièno zbrajati.

# Izraèun pomiènog zbroja CpG dinukleotida. 
# Prvo pomoæu vlastite funkcije only.CG sve CG dinukleotide
# zamjenjujem s 12, a onda sve 2 s 0. To radim da bi saèuvao
# poredak nukeotida (zato stavljam dvoznamenkasti broj), a zatim
# zamjenjujem 2 s 0 da bi mogao "prebrojati" CG zbrajanjem jedinica.
# Nakon funkcije je praktièki svaki CpG dinukleotid zamijenjen
# sa 1, ali tako da iza slijedi 0 da bi poredak nukleotida ostao
# saèuvan. Ostali nukleotidi su zamijenjeni s 0. 
# Zatim cijepam character vektor duljine 1 na vektor u kojem
# je svaki character odvojen. Na kraju na dobivenim brojevima
# funkcijom move.sum pomièno zbrajam CpG dinukleotide, a kako
# funkcija vraæa objekt klase Time-Series, sve pretvaram u brojeve
only.CG <- function(x){
  y <- gsub("CG", 12, x)
  y <- gsub(2, 0, y)
  y <- gsub("[GATC]", 0, y)
  return(y)
}
CpG.num.CG <- only.CG(CpG)
CpG.num.CG <- unlist(strsplit(CpG.num.CG, split = ""))
CpG.num.CG <- as.numeric(move.sum(CpG.num.CG))

# Izraèun pomiènog zbroja C nukleoida. Slièno kao i kod CpG
# dinukleotida, koristim isti princip i move.sum funkciju.
only.C <- function(x){
  y <- gsub("C", 1, x)
  y <- gsub("[GAT]", 0, y)
  return(y)
}
CpG.num.C <- only.C(CpG)
CpG.num.C <- unlist(strsplit(CpG.num.C, split = ""))
CpG.num.C <- as.numeric(move.sum(CpG.num.C))
  
# Izraèun pomiènog zbroja G nukleoida (sve isto kao iznad samo
# za G nukleotid)
only.G <- function(x){
  y <- gsub("G", 1, x)
  y <- gsub("[CAT]", 0, y)
  return(y)
}
CpG.num.G <- only.G(CpG)
CpG.num.G <- unlist(strsplit(CpG.num.G, split = ""))
CpG.num.G <- as.numeric(move.sum(CpG.num.G))


# Funkcija pomiènog prosjeka move.average:
# pomièe se po argumentu x za 1 uzimajuæi 100 èlanova i na 
# njima izraèunava prosjek. Koristi funkciju filter isto kao 
# move.sum (i iz istog razloga).
move.average <- function(x, n = 100){
  filter(x, rep(1/n, n), sides = 1) 
}

# Izraèun pomiènog prosjeka C+G nukleotida uz funkciju move.average.
# Princip je skoro isti kao i kod pomiènog zbrajanja, zamijenjujem
# sve G i C s 1, a sve A i T s O i pomoæu funkcije move.average
# raèunam pomièni prosjek. 
only.C.plus.G <- function(x){
  y <- gsub("[GC]", 1, x)
  y <- gsub("[AT]", 0, y)
  return(y)
}
CpG.avg.CG <- only.C.plus.G(CpG)
CpG.avg.CG <- unlist(strsplit(CpG.avg.CG, split = ""))
CpG.avg.CG <- as.numeric(move.average(CpG.avg.CG))


# tablica u kojoj su kolone po redu: 
# - pomièni zbroj CpG
# - pomièni zbroj C
# - pomièni zbroj G
# - izraèun pomiènog prosjeka Obs/Exp CpG
# - pomièni prosjek C+G nukleotida
CpG.df <- cbind(CpG.num.CG, CpG.num.C, CpG.num.G)
CpG.df <- data.frame(CpG.df)
CpG.df$"obs/exp" <- (CpG.df$CpG.num.CG/(CpG.df$CpG.num.C*CpG.df$CpG.num.G))*100
CpG.df <- cbind(CpG.df, CpG.avg.CG*100)
colnames(CpG.df) <- c("Number of CpG", "Number of C", 
                      "Number of G", "Obs/Exp", "%C+G")