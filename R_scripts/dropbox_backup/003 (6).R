CTCF <- read.table("K562-DS11247.peaks.fdr0.01.hg19.bed", stringsAsFactors = F, sep = "\t")
colnames(CTCF) <- c("chrom", "chromStart", "chromEnd", "name", "score",	"strand",
                    "signalValue", "pValue", "qValue", "peak") 

# Poredam data.frame po p-vrijednosti i uzimam 500 najveæih (500 prvih redova)

CTCF <- CTCF[with(CTCF, order(pValue, decreasing = T)), ] 
CTCF <- CTCF[1:500, ]


# Radim data.frame koji ima samo ime kromosoma i poèetne i krajnje pozicije 
# top 500 CTCF vezujuæih pozicija. Radim srednju vrijednost tih start i end 
# pozicija da dobijem centar sekvence. Onda tom centru oduzimam 200 da bi dobio
# start i dodajem 200 da bi dobio end sekvence. Te sekvence zapisujem u tablicu
# i sa UCSC table browsera skidam te sekvecnce. 

CTCF.position.data <- CTCF[, 1:3]
CTCF.position.data$center <- rowMeans(CTCF.position.data[2:3])
CTCF.position.data$chromStart <- CTCF.position.data$center - 200
CTCF.position.data$chromEnd <- CTCF.position.data$center + 200
CTCF.position.data <- CTCF.position.data[-4]

# Pomoæu Biostringov paketa i naredbe getSeq skidam sekvence po data.frameu 
# koji sam napravio iznad
library("Biostrings")
library("BSgenome.Hsapiens.UCSC.hg19")
library("seqLogo")

CTCF.data <- getSeq(Hsapiens, CTCF.position.data$chrom, CTCF.position.data$chromStart, 
                    CTCF.position.data$chromEnd)
CTCF.data <- as.character(CTCF.data)

# Prvo sam napisao manje-više sve funkcije koje koristim u Gibbsovom algoritmu:

# Uzimamo jednu po jednu sekvencu koju izdvojimo iz ostalih sekvenci
# Uzmemo nasumièan poèetak motiva za svaku od preostalih sekvenci
# i napravimo vektore motiva. Takoðer napravimo i vektor s nukleotidima 
# koji nisu u motivu. Sve skupa je vezano u data.frame i podijeljeno 
# u dvije funkcije 

RandomMotiffs <- function(sequence, motiff.length){
  N <- length(sequence)
  rand.motiff.start <- sample(nchar(sequence[1]) - motiff.length + 1, 
                              N, replace = T)
  sequence.df <- data.frame(sequence, rand.motiff.start)
  return(sequence.df)
}

MotiffsSequence <- function(sequence.df, motiff.length){
  rand.motiff <- mapply(substr, sequence.df[, 1],
                        sequence.df[, 2], sequence.df[, 2] + motiff.length - 1)
  sequence.df$motiff <- rand.motiff
  sequence.df$not.motiff <- sequence.df$sequence
  sequence.df$not.motiff <- mapply(sub, sequence.df[, 3], "", sequence.df[, 4])
  return(sequence.df)
}  


# Radimo position frequency matrix (PFM) - za svaku poziciju u motivu brojimo
# koliko ima razlièith baza. Prvo pojedine motive iz vektora pocijepam na
# pojedinaène baze. Zatim radim matricu koja ima u stupcima ima poziciju u motivu, 
# a u redovima pojedinaène motive. Zatim pomoæu funkcije nucl.column.freq
# brojim koliko u svakom stupcu imam pojedinih nukleotida. Primjenjujem 
# funkciju na stupce matrice i dobivam PFM.

PositionFrequencyMatrix <- function(sequence.df){
  
  rand.motiff <- strsplit(sequence.df$motiff, split = "")
  rand.motiff.df <- t(data.frame(rand.motiff))
  
  nucl.column.freq <- function(df){
    nucleotides <- c("A","C","G","T")
    vapply(nucleotides, function(x) length(grep(x, df)), integer(1))
  }
  
  PFM <- apply(rand.motiff.df, 2, nucl.column.freq)  
  
  # U 0 stupac PFM stavljam prebrojane nukelotide iz svih ne-motiva
  if(ncol(sequence.df) > 1){
    not.motiff <- unlist(strsplit(sequence.df$not.motiff, split = ""))
    PFM0 <- nucl.column.freq(not.motiff)
    PFM <- cbind(PFM0, PFM)
    colnames(PFM) <- as.character(c(0 : nchar(sequence.df$motiff[1])))  
  } else{
    colnames(PFM) <- as.character(c(1 : nchar(sequence.df$motiff[1])))  
  }
  
  return(PFM)
}


# Slijedeæe radim frequency count matrix po formuli:
# q(i, j) = ( c(i, j) + b(i) ) / (N - 1 + B)
# c(i, j) je broj u i-tom redu i j-tom stupcu PFM
# b(i) je konstanta za svaki nukleotid koju dodajemo kako bi dobili pseudocountove
# za vrijednosti gdje je broj nukleotida na toj poziciji 0 (kako nam kasnije
# vjerojatnost ne bi bila 0), uzeo sam da je on 0.5 za svaki nukleotid (našao
# na netu da je to okej vrijednost)
# N je broj sekvenci koje testiramo
# B je zbroj svih b(i), tj. u našem sluèaju 2

FrequencyCountMatrix <- function(PFM, N){
  
  freq.count <- function(df, N){
    q <- (df + 0.5)/(N + 2)  
  }
  
  FCM <- apply(PFM, c(1, 2), freq.count, N)  
  
  # U prvi stupac treba staviti podatke o background frequency countovima po formuli:
  # q (0, j) = ( c(0, j) + b(j) ) / ( sum(c(0, j)) + B )
  
  BFC <- PFM[, 1]
  BFC <- (BFC + 0.5) / (sum(BFC)+ 2)
  FCM[, 1] <- BFC
  return(FCM)  
}


# Tražim PWM za sekvencu koju sam izdvojio na poèetku
# Prvo napravim iteracije duljine motiva za tu sekvencu (funkcija ChosenSequenceMotiffIterations)

ChosenSequenceMotiffIterations <- function(chosen.seq, motiff.length = (ncol(my.FCM) - 1)){
  motiff.start <- 1 : (nchar(chosen.seq) - motiff.length + 1)
  chosen.seq.df <- data.frame(rep(chosen.seq, length(motiff.start)),
                              motiff.start)
  chosen.seq.iter <- mapply(substr, chosen.seq.df[, 1], chosen.seq.df[, 2], 
                            chosen.seq.df[, 2] + motiff.length - 1)
  return(chosen.seq.iter)
}


# Slijedeæa funkcija izraèunava PWM za pojedinaènu iteraciju dobivenu iznad.
# Funkcija uzima jednu iteraciju i raèuna weight prema FCM. Iteraciju podijelim
# na pojedinaène nukleotide i stavim ih u data.frame.
# Ideja je da taj nukleotid zamijenim brojem koji odgovara broju reda za 
# taj nukleotid u FCM (subset.GATC). Kad to napravim za sve nukleotide u df-u
# u jednoj varijabli imam koordinate redova u FCM-u u kojima se nalaze 
# nukleotidi za tu iteraciju. Koordinate stupaca su jednostavno brojevi od 
# 2 kraja FCM-a. To vežem u data.frame (chosen.seq.iter.coor). 
# Tada po tom data.frame-u s koordinatama vadim brojeve iz FCM. Da bi dobio
# position weight matrix za sve iteracije pokreæem funkciju sa sapply, 
# a rezultate spremam u data.frame. 

PWM.one.iteration <- function(chosen.seq.iter, FCM1 = my.FCM){
  
  mot.length <-  ncol(FCM1) - 1
  row.names(FCM1) <- c(1 : 4)
  chosen.seq.iter <- strsplit(chosen.seq.iter[1], split = "")
  chosen.seq.iter.df <- t(data.frame(chosen.seq.iter))
  rownames(chosen.seq.iter.df) <- NULL
  colnames(chosen.seq.iter.df) <- 1 : ncol(chosen.seq.iter.df)
  
  subset.GATC <- function(df){
    df <- gsub("A", 1, df)
    df <- gsub("C", 2, df)
    df <- gsub("G", 3, df)
    df <- gsub("T", 4, df)
  }
  
  chosen.seq.iter.rows <- apply(chosen.seq.iter.df, c(1, 2), subset.GATC)
  chosen.seq.iter.rows <- as.numeric(t(chosen.seq.iter.rows))
  chosen.seq.iter.coor <- rbind(chosen.seq.iter.rows, c(2 : (mot.length + 1)))
  
  prod1 <- 1
  prod2 <- 1
  
  for (i in 1 : mot.length){
    x <- chosen.seq.iter.coor[, i]
    prod1 <- prod(prod1, FCM1[x[1], x[2]])
    prod2 <- prod(prod2, FCM1[chosen.seq.iter.rows[i], 1])
  }
  div <- prod1 / prod2
  return(div)
}


# PWM je potrebno normalizirati. Svaka vrijednost u PWM se dijeli sa sumom
# svih vrijednosti.

NormalisePMW <- function(PWM){
  norm.PWM <- PWM / sum(PWM)
  return(norm.PWM)
}


# U PWM sada imamo probability weight za izabranu sekvencu. Prema tim vrijednostima
# biramo mjesto na kojem poèinje motiv u toj sekvenci.

MotiffStartChosenSeq <- function(PWM){
  motiff.start <- sample(ncol(PWM), 1, prob = PWM[1, ])
  return(motiff.start)
}


# Sam Gibbsov algoritam sam napravio ovako:
# Uzmem random sekvencu. Za ostale sekvence izaberem random motive. 
# Prema tim motivima napravim PFM i FCM. Prema FCM napravim weight matrix 
# za prije izabranu sekvencu. Prema tom weight matrixu izaberem motiv u izabranoj sekvenci. 
# Ponovno izaberem random sekvencu (drukèiju nego prošli put). 
# Za ostale sekvence izaberem random motive, ali za sekvencu koju sam prethodni 
# put izabrao stavim prošli put izabrani motiv. 
# Ponovno izraèunam PFM i FCM te weight za izabranu sekvencu. Prema weight 
# izaberem motiv za tu sekvencu. 
# Ponavljam postupak dok ne proðem sve sekvence puta. Napravim data.frame
# koji ima sekvence i mjesto poèetaka motiva koje sam dobio na taj naèin. 

# Zatim ponavljam postupak, ali umjesto nasumiènih poèetaka motiva koristim
# podatke za poèetak motiva dobivene u prvom koraku. To ponavljam n puta
# i svaki put motivi su sve precizniji (bliže zajednièkom motivu).

# Nažalost kad se radi na cijeloj sekvenci algoritam je dosta spor, pogotovo
# ako je n velik :(


motiff.length <- 14
n <- 10
my.full.sequence <- CTCF.data[1:5]
  
for (i in 1 : n){
  
  old.chosen.seq <- NULL
  all.chosen.seq <- NULL
  all.motiff.start <- NULL
  motiff.start <- NULL
  my.motiff.df <- NULL
  
  while(TRUE){
    
    while(TRUE){
      chosen.seq <- my.full.sequence[sample(length(my.full.sequence), 1)] 
      if(chosen.seq %in% all.chosen.seq == FALSE){
        break
      }
    }
    
    all.chosen.seq <- c(all.chosen.seq, chosen.seq)
    
    if(i == 1){
      my.sequence <- setdiff(my.full.sequence, chosen.seq)
      my.sequence.df <- RandomMotiffs(my.sequence, motiff.length)
    }
    
    if(i > 1){
      my.sequence.df <- my.sequence.df.1
      my.sequence.df <- my.sequence.df[my.sequence.df$sequence != chosen.seq, ]
    }
    
    my.sequence.df <- MotiffsSequence(my.sequence.df, motiff.length)
    my.sequence.df[my.sequence.df$sequence %in% old.chosen.seq, 2] <- motiff.start
    
    my.PFM <- PositionFrequencyMatrix(my.sequence.df)
    my.FCM <- FrequencyCountMatrix(my.PFM, length(my.full.sequence) - 1)
    
    my.chosen.seq.iter <- ChosenSequenceMotiffIterations(chosen.seq)
    my.PWM <- t(data.frame(sapply(my.chosen.seq.iter, PWM.one.iteration)))
    row.names(my.PWM) <- NULL
    my.PWM <- NormalisePMW(my.PWM)
    
    motiff.start <- MotiffStartChosenSeq(my.PWM)
    all.motiff.start <- c(all.motiff.start, motiff.start)
    if(length(all.motiff.start) == length(my.full.sequence)){
      break
    }
    
    old.chosen.seq <- chosen.seq
  }
  my.sequence.df.1 <- data.frame(all.chosen.seq, all.motiff.start)
  names(my.sequence.df.1) <- c("sequence", "rand.motiff.start")
}

consensus.motiffs <- my.sequence.df.1
consensus.motiffs <- MotiffsSequence(consensus.motiffs, motiff.length)
consensus.motiffs.PFM <- PositionFrequencyMatrix(consensus.motiffs)
consensus.motiffs.FCM <- FrequencyCountMatrix(my.PFM, length(my.full.sequence) - 1)
consensus.motiffs.PWM <- log(consensus.motiffs.FCM / 0.5)
consensus.motiffs.PWM <- apply(consensus.motiffs.PWM, 2, NormalisePMW)
consensus.motiffs.PWM <- consensus.motiffs.PWM[, 2 : motiff.length + 1]
seqLogo(consensus.motiffs.PWM)

# Nažalost nisam stigao napraviti ostatak zadatka (a i nisam baš siguran
# kako bi našao 3 najsnažnija motiva :) )