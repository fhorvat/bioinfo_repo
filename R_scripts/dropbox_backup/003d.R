CpG1 <- scan("141107_CpG", what = "character", sep = "\n")
x <- grep(">", CpG1)
CpG1 <- CpG1[-x]
CpG <- paste(CpG1,  sep="", collapse="")

only.CpG <- function(x){
  y <- gsub("CG", 12, x)
  y <- gsub(2, 0, y)
  y <- gsub("[GATC]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

only.C <- function(x){
  y <- gsub("C", 1, x)
  y <- gsub("[GAT]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

only.G <- function(x){
  y <- gsub("G", 1, x)
  y <- gsub("[CAT]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

only.CG <- function(x){
  y <- gsub("[GC]", 1, x)
  y <- gsub("[AT]", 0, y)
  y <- unlist(strsplit(y, split = ""))
  y <- as.numeric(y)
  return(y)
}

num.CpG <- only.CpG(CpG)
num.C <- only.C(CpG)
num.G <- only.G(CpG)
num.G.C <- only.CG(CpG)

n <- 247278
a <- 0
cnt <- rep(0, 100)

while (n < 247378){
  a <- a + 1
 
  x <- sum(num.CpG[(1 + n) : (200 + n)])
  y <- sum(num.C[(1 + n) : (200 + n)])
  z <- sum(num.G[(1 + n) : (200 + n)])
  
  Obs.Exp <- (x/(y*z))*200
  mov.avg.G.C <- sum(num.G.C[(1 + n) : (200 + n)]) / 200
  
  if (mov.avg.G.C > 50 & Obs.Exp > 0.6){
    cnt[a] <- n
  }
  
  n <- n + 1
}
cnt
