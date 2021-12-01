# CpG otoci su regije u DNA s visokim sadržajem 
# C+G nukleotida u kojima se èesto pojavljuju CpG
# dinukleotidi (s obzirom na ukupnu kolièinu DNA). 
# Inaèe je DNA kralješnjaka siromašna dinukleotidom
# CpG. 

# Definicija iz rada: CpG otok je DNA sekvenca za koju vrijedi:
# - prosjeèni broj G+C baza u sekvenci veæi od 0.5
# - omjer primjeæenih CpG/oèekivanih CpG dinukleotida
# (Obs/Exp CpG) je veæi od 0.6
# - sekvenca je dugaèka najmanje 200 bp

# Obs/Exp CpG = ((number of CpG)/(number of C * number of G)) * N
# gdje je N ukupan broj nukleotida u analiziranoj sekvenci

# U radu promatraju sekvencu pomoæu "prozora" od 200 bp (N=200)
# koji se pomièe jednu po jednu bazu pa sam i ja tako radio. 

CpG <- scan("141107_CpG", what = "character", sep = "\n")
x <- grep(">", CpG)
CpG <- CpG[ - x]
CpG <- paste(CpG,  sep = "", collapse = "")

# Funkcija koja kotrljajuæe zbraja pomoæu funkcije filter koja
# inaèe radi na time.series podacima, ali mi je korisna jer pomoæu
# nje mogu pomièno zbrajati s prozorom duljine n. Odluèio
# sam korisiti funkciju filter nakon što su mi se petlje na
# ovoliko dugaèkim podacima vrtile nekoliko desetaka minuta. 
move.sum <- function(x, n = 200){
  y <- as.numeric(filter(x, rep(1, n), sides = 1))
  return(y)
}

# Da bi mogao kotrljajuæe zbrajati nukleotide pomoæu funkcije
# move.sum, moram ih pretvoriti u brojeve. To radim pomoæu 
# slijedeæe 4 funkcije:

# Funkcija koja pretvara dinukleotid CG u broj 1. Prvo u stringu
# tražim "CG" i zamjenjujem ga brojem 10. Koristim dvoznamenkasti
# broj zato da bi saèuvao poredak nukleotida (kad bih stavio 
# samo 1 narušio bi poredak). Zatim 2 sve preostale baze pretvaram
# u 0. 
only.CpG <- function(x){
  y <- gsub("CG", 10, x)
  y <- gsub("[GATC]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

# Slièno kao i funkcija iznad, jednostavno sve C zamijenjujem
# s 1, a ostale nukleotide s 0
only.C <- function(x){
  y <- gsub("C", 1, x)
  y <- gsub("[GAT]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

# Isto kao i funkcija iznad samo za G nukleotid
only.G <- function(x){
  y <- gsub("G", 1, x)
  y <- gsub("[CAT]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

# Slièno kao i ostale funkcije, G i C zamijenjujem s 1, A i T
# s 0.
only.GC <- function(x){
  y <- gsub("[GC]", 1, x)
  y <- gsub("[AT]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

# Napravim 4 varijable s 4 gornje funkcije, svaka varijabla
# ima drugu kombinaciju nukleotida zamijenjenu brojevima. 
a <- only.CpG(CpG)
b <- only.C(CpG)
d <- only.G(CpG)
e <- only.GC(CpG)

# Pomoæu funkcije move.sum pomièno zbrajam pojedine nukelotide
# ili dinukleotide. Dobijem varijablu u kojoj su mi zbrojeni
# nukleotidi s prozorom 200 i pomakom od 1. 
a <- move.sum(a)
b <- move.sum(b)
d <- move.sum(d)
e <- move.sum(e)

# Varijablu e u kojoj su mi zbojeni G i C nukleotidi dijelim sa 
# 200 kako bi dobio prosjeènu vrijednost u G i C nukleotida
# u sekvenci od 200 bp. 
# Izraèunavam omjer primjeæenih CpG i oèekivanih CpG dinukleotida
# (Obs/Exp CpG)
e <- e/200
Obs.Exp <- (a / (b * d)) * 200

# tablica: prosjeèan broj C+G nukleotida i Obs/Exp CpG     
tbl <- data.frame(cbind(e, Obs.Exp))

# Da bi našao CpG otoke prvo izdvajam sve toèke koje ispunjavaju:
# - Obs/Exp CpG > 0.6
# - prosjeèni broj G+C baza (e) u sekvenci veæi od 0.5
tbl <- subset(tbl, Obs.Exp > 0.6 & e > 0.5)

# Dodajem novi stupac u kojem su pozicije toèaka u kromosomu.
# Broju reda dodajem poèetnu poziciju sekvence. Ta poèetna pozicija
# nije jednaka broju koji je zadan u zadataku, naime kada na UCSC 
# table browseru stavim poèetak sekvence koji je zadan u zadatku
# UCSC mi vraæa  sekvencu koja poèinje oko 80 000 bp prije. Još 
# oduzimam 200 koliko je dugaèak prozor. 

tbl$position <- as.numeric(row.names(tbl)) + 22742858 - 200

# Pretpostavljam da sve toèke koje slijede jedna iza druge, 
# a zadovoljavaju gornja dva uvjeta pripadaju istom potencijalnom
# CpG otoku. Ne znam je li to toèno, ali mi je logièno pa 
# sam dalje tako radio. Sad želim dobiti pozicije na kromosomu
# svakog poèetka i kraja pojedinog otoka. 

# Raèunam razmak izmeðu slijedeæe i prethodne toèke (2 - 1, 
# 3 - 2, 4 - 3 itd.) i dobivene podatke stavljam u novu kolonu 
# u tablici. Tamo gdje je razmak izmeðu slijedeæe i prethodne
# toèke veæi od 1 poèinje novi otok.
tbl$position.2 <- c((tbl[2:nrow(tbl), 3] - tbl[1:(nrow(tbl)-1), 3]), 1)

# Radim novu varijablu u koju stavljam pozicije na kromosomu 
# gdje poèinju CpG otoci. 
start <- tbl[1, 3] 
start <- c(start, (tbl[tbl[ , 4] > 1, 4] + tbl[tbl[ , 4] > 1, 3]))

# Isto tako radim varijablu u kojoj su pozicije krajeva CpG otoka
# na kromosomu 
end <- tbl[tbl[ ,4] > 1, 3]
end <- c(end, tbl[nrow(tbl), 3])

# Stavljam pozicije poèetaka i krajeva otoka u novu tablicu.
# U novom stupcu raèunam duljinu CpG otoka
mydata <- data.frame(cbind(start, end))
mydata$length <- mydata$end - mydata$start

# Kako su CpG otoci po definiciji iz rada duži od 200 bp, 
# odbabujem sve sekvence koje su kraæe od 200 bp.
mydata <- subset(mydata, length > 200)

# Uèitavam podatke za CpG otoke sa UCSC: 
ucsc.data <- read.table("141107_CpG_islands", header = F, sep = "\t")
colnames(ucsc.data) <- c("bin",  "chrom",  "start",	"end",	
                         "name",	"length",	"cpgNum",	"gcNum",	
                         "perCpg",	"perGc",	"obsExp")
ucsc.data <- ucsc.data[ , c(-1, -2, -5, -7:-11)]

# Pomoæu podataka o poèetnom i krajnjem mjestu svakog otoka
# radim sekvencu brojeva za svaki od otoka (od start do end
# za svaki otok) i stavljam sve te brojeve u jedan vektor
mydata.seq <- mapply(":", mydata$start, mydata$end)
mydata.seq <- unlist(mydata.seq)

# Isto radim za podatke s UCSC-a
ucsc.data.seq <- mapply(":", ucsc.data$start, ucsc.data$end)
ucsc.data.seq <- unlist(ucsc.data.seq)

# Gledam koje pozicije u kromosomu se pojavljuju u obje sekvence i
# pomoæu which i %in% i stavljam ih u varijablu overlap
overlap <- ucsc.data.seq[which(ucsc.data.seq %in% mydata.seq)]

# Da bi izraèunao poèetak i kraj tih sekvenci prvo ih stavljam
# u data.frame pa u novu kolonu stavljam razliku slijedeæe pozicije
# u kromosomu od prethodne. 
overlap.df <- data.frame(overlap)
overlap.df$razlika <- c((overlap.df[2:nrow(overlap.df), 1] 
                         - overlap.df[1:(nrow(overlap.df)-1), 1]), 0)

# Pomoæu podataka iz gornjeg data.frame-a u varijable start i end
# stavljam podatke o poèetku i kraju svake sekvence.
start <- overlap.df[1, 1] 
start <- c(start, (overlap.df[overlap.df[ , 2] > 1, 2] + overlap.df[overlap.df[ , 2] > 1, 1]))
end <- overlap.df[overlap.df[ , 2] > 1, 1]
end <- c(end, overlap.df[nrow(overlap.df), 1])

# U novi data.frame stavljam podatke o poèetnim i krajnjim
# pozicijama preklapanja CpG otoka u kromosomu
CpG.overlap <- data.frame(cbind(start, end))
CpG.overlap$length <- CpG.overlap$end - CpG.overlap$start

# Mojih 97 CpG otoka se preklapa s 28 CpG otoka dobivenih s UCSC-a
# u jedanaest sekvenci koje su prikazane u data.frame-u:
CpG.overlap

# sve sekvence poredane po poziciji u kromosomu su prikazane u 
# data.frameu all.data
CpG.overlap$id <- "overlap"
mydata$id <- "mydata"
ucsc.data$id <- "ucsc"

all.data <- rbind(CpG.overlap, ucsc.data, mydata)
all.data <- all.data[-3]
all.data <- all.data[order(all.data[, 1], all.data[, 2]), ]

# Confusion matrix pozicija u kromosomu
confusion.matrix <- function(x, y){
  both <-  union(x, y)
  actual <-  both %in% x
  predicted <-  both %in% y
  return(table(actual, predicted))
}
conf.mat <- confusion.matrix(ucsc.data.seq, mydata.seq)
conf.mat

# Data.frame s podacima iz confusion matrixa
conf.mat.df <- data.frame(conf.mat)
conf.mat.df <- data.frame(lapply(conf.mat.df, as.character), 
                          stringsAsFactors = FALSE)
conf.mat.df$description <- c("TN", "FN", "FP", "TP") 

# Raèunanje: 
# Sensitivity  = TP / (TP + FN)
# Specificity = TN / (FP + TN)
# False discovery rate = FP / (FP + TP)
# Accuracy = (TP + TN) / (TP + FN + FP + TN)
#   TP = true positive   TN = true negative
#   FP = false positive  FN = false negative 

TN <- as.numeric(conf.mat.df[1, 3])
FN <- as.numeric(conf.mat.df[2, 3])
FP <- as.numeric(conf.mat.df[3, 3])
TP <- as.numeric(conf.mat.df[4, 3])

# Sensitivity
sens  = TP / (TP + FN)
sens

# Specificity
spec = TN / (FP + TN)
spec

# False discovery rate
fdr = FP / (FP + TP)
fdr

# Accuracy
acc = (TP + TN) / (TP + FN + FP + TN)
acc
