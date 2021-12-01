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

FrequencyCountMatrix <- function(PFM, N = length(my.full.sequence), BFC.t.f){
  
  freq.count <- function(df, N){
    q <- (df + 0.5)/(N - 1 + 2)  
  }
  
  FCM <- apply(PFM, c(1, 2), freq.count, N)  
  
  # U prvi stupac treba staviti podatke o background frequency countovima po formuli:
  # q (0, j) = ( c(0, j) + b(j) ) / ( sum(c(0, j)) + B )
  
  if (BFC.t.f == T){
    BFC <- PFM[, 1]
    BFC <- (BFC + 0.5) / (sum(BFC)+ 2)
    FCM[, 1] <- BFC
    return(FCM)  
  }
  
  if (BFC.t.f == F){
    return(FCM)  
  }
  
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

SplitMotiff <- function(chosen.seq, motiff.startm, motiff.length){
  N <- length(chosen.seq)
  sequence.df <- data.frame(chosen.seq, motiff.start)
  motiff <- mapply(substr, sequence.df[, 1],
                   sequence.df[, 2], sequence.df[, 2] + motiff.length - 1)
  return(motiff)
}

# Uzmem random sekvencu. Za ostale sekvence izaberem random motive. 
# Prema tim motivima napravim PFM i FCM. Prema FCM napravim weight matrix 
# za prije izabranu sekvencu. 
# Prema tom weight matrixu izaberem motiv u izabranoj sekvenci i spremim ga. Ponovno
# Ponovno izaberem random sekvencu (drukèiju nego prošli put). 
# Za ostale sekvence izaberem random motive, ali za sekvencu koju sam prethodni 
# put izabrao stavim prošli put izabrani motiv. 
# Ponovno izraèunam PFM i FCM te weight za izabranu sekvencu. Prema weight 
# izaberem motiv za tu sekvencu. Spremim taj motiv. 
# Ponavljam postupak i puta. Spremam sve izabrane motive u data.frame
# Iz tog data.frame-a radim consensus frequency count matrix. 
  

motiff.length <- 7
my.full.sequence <- c("ACCATGACAG", "GAGTATACCT", "CATGCTTACT", "CGGAATGCAT")
old.chosen.seq <- NULL
my.motiff.df <- NULL
motiff.start <- NULL

for (i in 1:500){
  
  while(TRUE){
    chosen.seq <- my.full.sequence[sample(length(my.full.sequence), 1)] 
    if(chosen.seq %in% old.chosen.seq == FALSE){
      break
    }
  }
  
  my.sequence <- setdiff(my.full.sequence, chosen.seq)
  my.sequence.df <- RandomMotiffs(my.sequence, motiff.length)
  my.sequence.df <- MotiffsSequence(my.sequence.df, motiff.length)
  my.sequence.df[my.sequence.df$sequence %in% old.chosen.seq, 2] <- motiff.start
  
  my.PFM <- PositionFrequencyMatrix(my.sequence.df)
  my.FCM <- FrequencyCountMatrix(my.PFM, BFC.t.f = T)
  
  my.chosen.seq.iter <- ChosenSequenceMotiffIterations(chosen.seq)
  my.PWM <- t(data.frame(sapply(my.chosen.seq.iter, PWM.one.iteration)))
  row.names(my.PWM) <- NULL
  my.PWM <- NormalisePMW(my.PWM)
  
  motiff.start <- MotiffStartChosenSeq(my.PWM)
  my.motiff <- SplitMotiff(chosen.seq, motiff.start, motiff.length)
  
  if (is.null(my.motiff.df) == TRUE){
    my.motiff.df <- t(data.frame(my.motiff))  
  } else{
    my.motiff.df <- rbind(my.motiff.df, my.motiff)
  }
  
  old.chosen.seq <- chosen.seq
}

my.motiff.df <- data.frame(my.motiff.df, stringsAsFactors = F, row.names = NULL)
names(my.motiff.df) <- "motiff"

my.motiff.PFM <- PositionFrequencyMatrix(my.motiff.df)
my.motiff.FCM <- FrequencyCountMatrix(my.motiff.PFM, BFC.t.f = F)
my.motiff.FCM
PWM(my.motiff.PFM)
?seqLogo
