setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/Projekti/test/adventOfCode")

input1 <- read.table("day03_input_1.txt", sep = "", stringsAsFactors = F)

possibleTriangles <- function(input_row){
  input_row <- 
  input_comb_sum <- rowSums(t(combn(input_row, 2)))
  possible <- all(unlist(lapply(X = input_row, FUN = function(X) input_comb_sum > X)))
  return(possible)
}

sum(apply(input1, 1, possibleTriangles))

# part 2
input2 <- matrix(unname(unlist(input1)), ncol = 3, byrow = T)
sum(apply(input2, 1, possibleTriangles))

# solution 2
sum(apply(read.table("day03_input_1.txt", sep = ""), 1, function(x) sum(x) > 2 * max(x)))
sum(apply(matrix(unlist(read.table("day03_input_1.txt", sep = "")), ncol = 3, byrow = T), 1, function(x) sum(x) > 2 * max(x)))
