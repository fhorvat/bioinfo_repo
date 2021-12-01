MRSA <- scan("MRSA252.dna", what = "character", sep = "\n")
x <- grep(">", MRSA)
MRSA <- MRSA[-x]

# Funkcija koja uvodi netoène baze u èitanja. Stavio sam da A zamjenjuje s T i obrnuto, 
# a G s C i obrnuto. Da ne bi za svaku sekvencu morao vrtiti petlju onoliko puta koliko
# je sekvenca duga i za svaku bazu raditi sample (1000, 3), napravio sam sample
# koji izvuèe onoliko brojeva koliko je duga sekvenca iz brojeva od 1 do 1000
# i onda gledam je li izvukao broj 3 ili manji i koliko takvih brojeva je izvukao.
# Onda vrtim petlju od 1 do broja brojeva manjih ili jednak od 3 i mijenjam toliko 
# baza. 
       
incorrect.read.fun <- function(sequence){
  incorrect.read <- sample(1000, nchar(sequence), replace = T)
  incorrect.pos <- grep("^[1-3]$", incorrect.read)
  if (length(incorrect.pos) > 0){
    for (l in 1 : length(incorrect.pos)){
      if (substr(sequence, incorrect.pos[l], incorrect.pos[l]) == "a"){
        substr(sequence, incorrect.pos[l], incorrect.pos[l]) <- "t"
      } 
      if (substr(sequence, incorrect.pos[l], incorrect.pos[l]) == "t"){
        substr(sequence, incorrect.pos[l], incorrect.pos[l]) <- "a"
      }
      if (substr(sequence, incorrect.pos[l], incorrect.pos[l]) == "c"){
        substr(sequence, incorrect.pos[l], incorrect.pos[l]) <- "g"
      }
      if (substr(sequence, incorrect.pos[l], incorrect.pos[l]) == "g"){
        substr(sequence, incorrect.pos[l], incorrect.pos[l]) <- "c"
      }  
    }
  }
  return(sequence)    
}

my.kmer.coverage <- function(sequence, kLen, coverage){
  
  kmers.lst <- NULL
  kmers.lst.full <- NULL
  
  for (i in 1 : ((coverage * nchar(sequence))/100)){
    seq.length <- round(rnorm(1, mean = 100, sd = 6))
    start.pos <- sample(nchar(sequence), 1)
    all.kmers <- substr(sequence, start.pos, (start.pos + seq.length))
    all.kmers.incorrect <- incorrect.read.fun(all.kmers)
    
    for (j in 1 : (nchar(all.kmers.incorrect) - kLen)){
      kmers.lst <- c(kmers.lst, substr(all.kmers.incorrect, 0 + j, kLen + j))  
    }
    
    kmers.lst.full <- c(kmers.lst.full, kmers.lst)
    kmers.lst <- NULL
  }
  
  kmers.lst.full.df <- as.data.frame(table(table(kmers.lst.full)), 
                                     stringsAsFactors = F)
  return(plot(x = kmers.lst.full.df[, 1], y = kmers.lst.full.df[, 2], pch = 20))
}

MRSA1 <- MRSA[1:700]
MRSA1 <- unlist(strsplit(MRSA1, split = ""))
MRSA1 <- paste(MRSA1,  sep = "", collapse = "")
my.kmer.coverage(sequence = MRSA1, kLen = 20, coverage = 15)

# Svi kmerovi koji se pojavljuju jednom imaju grešku pri èitanju. 
# Kako je mutacija prilièno rijetka, a sekvenca velika i dobro prekrivena 
# tijekom sekvenciranja (veliki coverage), mala je šansa 
# da se ista mutacija ponovi toèno na istom mjestu pa se zato isti kmerovi 
# s mutacijom ne ponavljaju. Na grafu su to kmerovi s vrlo malom vrijednosti 
# na x-osi (pojavljuju se jednom ili najviše par puta tijekom sekvenciranja) i 
# s vrlo visokom vrijednosti na  y-osi (ima puno takvih kmerova koji se ponavljaju
# mali broj puta).