library(stringr)
library(ggplot2)

setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/Projekti/test/adventOfCode")
input <- scan("day02_input_1.txt", what = "character")
input <- strsplit(input, split = "")
  
map <- setNames(c(1, -1, -1, 1), c("R", "L", "U", "D"))
keyboard <- matrix(1:9, nrow = 3, ncol = 3, byrow = T)

keyboardMove <- function(X){
  
  move <- map[input_sub[X]]
  
  if(names(move) %in% c("L", "R") & position[1] >= 1 & position[1] <= 3){
    position_holder <- position[1] + unname(move)
    if(position_holder != 0 & position_holder != 4){
      position[1] <<- position[1] + unname(move)
    }
  }
  
  if (names(move) %in% c("D", "U") & position[2] >= 1 & position[2] <= 3){
    position_holder <- position[2] + unname(move)
    if(position_holder != 0 & position_holder != 4){
      position[2] <<- position[2] + unname(move)
    }
  }
}

position <- c(2, 2)
pin <- NULL
for(i in 1:length(input)){
  input_sub <- input[[i]]
  invisible(lapply(X = 1:length(input_sub), FUN = keyboardMove))
  pin <- c(pin, keyboard[position[2], position[1]])
}
pin 

# part 2
keyboard2 <- rbind(c(NA, NA, 1, NA, NA), 
                   c(NA, 2, 3, 4, NA), 
                   c(5, 6, 7, 8, 9), 
                   c(NA, "A", "B", "C", NA), 
                   c(NA, NA, "D", NA, NA))

keyboardMove2 <- function(X){
  
  move <- map[input_sub[X]]
  
  if(names(move) %in% c("L", "R") & position[1] >= 1 & position[1] <= 5){
    position_holder <- position[1] + unname(move)
    if(position_holder != 0 & position_holder != 6){
      if(!is.na(keyboard2[position[2], position_holder])){
        position[1] <<- position[1] + unname(move)
      }
    } 
  }

  
  if(names(move) %in% c("D", "U") & position[2] >= 1 & position[2] <= 5){
    position_holder <- position[2] + unname(move)
    if(position_holder != 0 & position_holder != 6){
      if(!is.na(keyboard2[position_holder, position[1]])){
        position[2] <<- position[2] + unname(move)
      }
    } 
  }
}

position <- c(1, 3)
pin <- NULL
for(i in 1:length(input)){
  input_sub <- input[[i]]
  invisible(lapply(X = 1:length(input_sub), FUN = keyboardMove2))
  pin <- c(pin, keyboard2[position[2], position[1]])
}
pin 





