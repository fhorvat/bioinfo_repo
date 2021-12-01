setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/Projekti/test/adventOfCode")

input1 <- read.fwf("day06_input_1.txt", widths = rep(1, 8), stringsAsFactors = F)
input_count <- apply(input1, 2, table)
paste0(rownames(input_count[apply(input_count, 2, which.max), ]), collapse = "")

# ugly one-liner
paste0(rownames(input_count[apply(apply(read.fwf("day06_input_1.txt", widths = rep(1, 8), stringsAsFactors = F), 2, table), 2, which.max), ]), collapse = "")

# part 2
paste0(rownames(input_count[apply(input_count, 2, which.min), ]), collapse = "")

# ugly one-liner part 2
paste0(rownames(input_count[apply(apply(read.fwf("day06_input_1.txt", widths = rep(1, 8), stringsAsFactors = F), 2, table), 2, which.min), ]), collapse = "")
