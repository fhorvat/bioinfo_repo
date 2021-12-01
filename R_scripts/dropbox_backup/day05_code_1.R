setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/Projekti/test/adventOfCode")
library(digest)

input <- "ojvtpuvg"
pass <- NULL
n <- 0
repeat{
  n <- n + 1
  md5_hash <- digest(paste0(input, n), algo = 'md5', serialize = FALSE)
  if(grepl("^00000", md5_hash)){
    pass <- c(pass, strsplit(md5_hash, split = "")[[1]][6])
  }
  if(length(pass) == 8){
    break
  }
}

paste0(pass, collapse = "")

# part 2
pass2 <- rep("*", 8)
n <- 0
count <- 0
repeat{
  n <- n + 1
  md5_hash <- digest(paste0(input, n), algo = 'md5', serialize = FALSE)
  if(grepl("^00000", md5_hash)){
    position <- strsplit(md5_hash, split = "")[[1]][6]
    if(!is.na(suppressWarnings(as.integer(position)))){
      if(as.integer(position) >= 0 & as.integer(position) <= 7){
        position <- as.integer(position) + 1
        if(pass2[position] == "*"){
          pass2[position] <- strsplit(md5_hash, split = "")[[1]][7]
          count <- count + 1
        }
      }
    }
}
  if(count == 8){
    break
  }
}
paste0(pass2, collapse = "")
