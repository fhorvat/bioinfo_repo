MRSA <- scan("MRSA252.dna", what = "character", sep = "\n")
x <- grep(">", MRSA)
MRSA <- MRSA[-x]

MRSA1 <- MRSA[1:60]
MRSA1 <- unlist(strsplit(MRSA1, split = ""))
MRSA1 <- paste(MRSA1,  sep = "", collapse = "")

base.cov <- 15
kLen <- 20
kmers.lst <- list()
kmers.lst1 <- list()
kmers.lst2 <- list()
inc.read.lst <- list()
kmers.lst.full1 <- list()
kmers.lst.full2 <- list()


for (i in 1 : (base.cov * nchar(MRSA1))){
  
  seq.length <- round(rnorm(1, mean = 100, sd = 6))
  start.pos <- sample(nchar(MRSA1), 1)
  sequence <- substr(MRSA1, start.pos, (start.pos + seq.length))
  
  incorrect.read <- sample(1000/seq.length, 1)
  
  if (incorrect.read <= 3){
    sequence1 <- sequence
    incorrect.read.pos <- sample(nchar(sequence1), 1)
    incorrect.read.seq <- substr(sequence1, incorrect.read.pos, incorrect.read.pos)
    
    if (incorrect.read.seq == "a"){
      substr(sequence1, incorrect.read.pos, incorrect.read.pos) <- "t"
    }
    if (incorrect.read.seq == "t"){
      substr(sequence1, incorrect.read.pos, incorrect.read.pos) <- "a"
    }
    if (incorrect.read.seq == "c"){
      substr(sequence1, incorrect.read.pos, incorrect.read.pos) <- "g"
    }
    if (incorrect.read.seq == "g"){
      substr(sequence1, incorrect.read.pos, incorrect.read.pos) <- "c"
    }
    
    for (j in 0 : (nchar(sequence1) - kLen)){
      kmers.lst1[j] <- list(substr(sequence1, 0 + j, kLen + j))  
    }
    kmers.lst.full1[i] <- list(kmers.lst1)
  }
  
  else{
    for (j in 0 : (nchar(sequence) - kLen)){
      kmers.lst2[j] <- list(substr(sequence, 0 + j, kLen + j))  
    }
    kmers.lst.full2[i] <- list(kmers.lst2)
  }
  
  kmers.lst1 <- list()
  kmers.lst2 <- list()
}

kmers.lst.full1 <- unlist(kmers.lst.full1)
kmers.lst.full2 <- unlist(kmers.lst.full2)

kmers.lst1.tbl <- as.data.frame(table(table(kmers.lst.full1)), stringsAsFactors = F)
kmers.lst2.tbl <- as.data.frame(table(table(kmers.lst.full2)), stringsAsFactors = F)
plot(x = kmers.lst1.tbl[, 1], y = kmers.lst1.tbl[, 2], pch = 20)
points(x = kmers.lst2.tbl[, 1], y = kmers.lst2.tbl[, 2], pch = 20, col = "red")
