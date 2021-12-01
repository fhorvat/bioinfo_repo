sequence <- c("ACCATGACAG", "GAGTATACCT", "CATGCTTACT",
              "CGGAATGCAT")
motiff <- "GATTACA"
motiff.length <- 7
N <- length(sequence) 

# Uzimamo jednu nasumiènu sekvencu koju izdvojimo iz ostalih sekvenci

chosen.seq <- sequence[sample(length(sequence), 1)]
sequence <- setdiff(sequence, chosen.seq)


# Uzmemo nasumièan poèetak motiva za svaku od preostalih sekvenci
# i napravimo vektor s motivima. Takoðer napravimo i vektor s nukleotidima 
# koji nisu u motivu
rand.motiff.start <- sample(nchar(sequence[1]) - motiff.length + 1, (N - 1), 
                            replace = T)
sequence.df <- data.frame(sequence, rand.motiff.start)
rand.motiff <- mapply(substr, sequence.df[, 1], sequence.df[, 2], 
                      sequence.df[, 2] + motiff.length - 1)
sequence.df$motiff <- rand.motiff
sequence.df$not.motiff <- sequence.df$sequence
sequence.df$not.motiff <- mapply(sub, sequence.df[, 3], "", sequence.df[, 4])


# Radimo position frequency matrix (PFM) - za svaku poziciju u motivu brojimo
# koliko ima razlièith baza. Prvo pojedine motive iz vektora pocijepam na
# pojedinaène baze. Zatim radim matricu koja ima u stupcima ima poziciju u motivu, 
# a u redovima pojedinaène motive. Zatim pomoæu funkcije nucl.column.freq
# brojim koliko u svakom stupcu imam pojedinih nukleotida. Primjenjujem 
# funkciju na stupce matrice i dobivam PFM.

rand.motiff <- strsplit(rand.motiff, split = "")
rand.motiff.df <- t(data.frame(rand.motiff))

nucl.column.freq <- function(df){
  nucleotides <- c("A","C","G","T")
  vapply(nucleotides, function(x) length(grep(x, df)), integer(1))
}

PFM <- apply(rand.motiff.df, 2, nucl.column.freq)


# U 0 stupac PFM stavljam prebrojane nukelotide iz svih ne-motiva

not.motiff <- unlist(strsplit(sequence.df$not.motiff, split = ""))
PFM0 <- nucl.column.freq(not.motiff)
PFM <- cbind(PFM0, PFM)
colnames(PFM) <- as.character(c(0 : motiff.length))


# Slijedeæe radim frequency count matrix po formuli:
# q(i, j) = ( c(i, j) + b(i) ) / (N - 1 + B)
# c(i, j) je broj u i-tom redu i j-tom stupcu PFM
# b(i) je konstanta za svaki nukleotid koju dodajemo kako bi dobili pseudocountove
# za vrijednosti gdje je broj nukleotida na toj poziciji 0 (kako nam kasnije
# vjerojatnost ne bi bila 0), uzeo sam da je on 0.5 za svaki nukleotid (našao
# na netu da je to okej vrijednost)
# N je broj sekvenci koje testiramo
# B je zbroj svih b(i), tj. u našem sluèaju 2

freq.count <- function(df, N){
  q <- (df + 0.5)/(N - 1 + 2)  
}

FCM <- apply(PFM, c(1, 2), freq.count, N)


# U prvi stupac treba staviti podatke o background frequency countovima po formuli:
# q (0, j) = ( c(0, j) + b(j) ) / ( sum(c(0, j)) + B )

BFC <- PFM[, 1]
BFC <- (BFC + 0.5) / (sum(BFC)+ 2)
FCM[, 1] <- BFC


# Tražim PWM za sekvencu koju sam izdvojio na poèetku
# Prvo napravim iteracije duljine motiva za tu sekvencu 
motiff.start <- 1 : (nchar(chosen.seq) - motiff.length + 1)
chosen.seq.df <- data.frame(rep(chosen.seq, length(motiff.start)),
                            motiff.start)
chosen.seq.iter <- mapply(substr, chosen.seq.df[, 1], chosen.seq.df[, 2], 
                          chosen.seq.df[, 2] + motiff.length - 1)

# Funkcija uzima jednu iteraciju i raèuna weight prema FCM za ostale sekvence
# Na kraju imam weight za sve iteracije. 

PWM.iterations <- function(chosen.seq.iter, mot.length = motiff.length, 
                           FCM1 = FCM){
  
  row.names(FCM1) <- c(1 : 4)
  chosen.seq.iter <- strsplit(chosen.seq.iter, split = "")
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
  chosen.seq.iter.coor.1 <- rbind(chosen.seq.iter.rows, c(2 : (mot.length + 1)))
  
  prod1 <- 1
  prod2 <- 1
  
  for (i in 1 : mot.length){
    x <- chosen.seq.iter.coor.1[, i]
    prod1 <- prod(prod1, FCM1[x[1], x[2]])
    prod2 <- prod(prod2, FCM1[chosen.seq.iter.rows[i], 1])
  }
  div <- prod1 / prod2
  return(div)
}

a <- t(data.frame(sapply(chosen.seq.iter, PWM.iterations)))

