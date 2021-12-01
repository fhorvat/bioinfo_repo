library(stringr)
library(ggplot2)

setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/Projekti/test/adventOfCode")
input <- scan("day01_input_1.txt", what = "character", sep = ",", strip.white = T)

# input <- c("R8", "R4", "R4", "R8")
sides_df <- data.frame(dir = c("N", "S", "E", "W"), 
                       L = c("W", "E", "N", "S"), 
                       R = c("E", "W", "S", "N"), 
                       stringsAsFactors = F)

dir <- "N"
input_world <- NULL
for(step in input){
  step_dir <- gsub("[[:digit:]]", "", step)
  step_n <- gsub("[[:alpha:]]", "", step)
  dir <- sides_df[sides_df$"dir" == dir, step_dir]
  input_world <- c(input_world, paste0(dir, step_n))
}

east_west <- input_world[grepl("E|W", input_world)]
east_west <- gsub("W", "-", east_west)
east_west <- gsub("E", "", east_west)
east_west <- as.numeric(east_west)
east_west_sum <- sum(east_west)

north_south <- input_world[grepl("S|N", input_world)]
north_south <- gsub("S", "-", north_south)
north_south <- gsub("N", "", north_south)
north_south <- as.numeric(north_south)
north_south_sum <- sum(north_south)

sum(abs(east_west_sum), abs(north_south_sum))

# part 2
input_world_coordinates <- gsub("W|S", "-", input_world)
input_world_coordinates <- gsub("E|N", "", input_world_coordinates)

input_world_df <- data.frame(input_world, 
                             coordinates = input_world_coordinates, 
                             stringsAsFactors = F)
x_coordinate <- NULL
y_coordinate <- NULL
for(i in 1:nrow(input_world_df)){
  
  if(grepl("E|W", input_world_df$input_world[i])){
    x_coordinate <- c(x_coordinate, input_world_df$coordinates[i])
  }else{
    x_coordinate <- c(x_coordinate, 0)
  }
  
  if(grepl("S|N", input_world_df$input_world[i])){
    y_coordinate <- c(y_coordinate, input_world_df$coordinates[i])
  }else{
    y_coordinate <- c(y_coordinate, 0)
  }

}

x_coordinate_sum <- cumsum(x_coordinate)
x_coordinate_sum <- x_coordinate_sum[seq(1, length(x_coordinate_sum), 2)]
x_coordinate_sum <- c(0, x_coordinate_sum)

y_coordinate_sum <- cumsum(y_coordinate)
y_coordinate_sum <- y_coordinate_sum[seq(2, length(y_coordinate_sum), 2)]
y_coordinate_sum <- c(0, 0, y_coordinate_sum)

coordinate_df <- data.frame(x_coordinate_sum, y_coordinate_sum)

step_1_x <- step_1_y <- step_2_x <- step_2_y <- step_x <- step_y <- NULL
for (i in 1:(length(x_coordinate_sum) - 2)){
  step_1_x <- x_coordinate_sum[i] : x_coordinate_sum[i + 1]
  step_1_y <- rep(y_coordinate_sum[i + 1], length(x_coordinate_sum[i]:x_coordinate_sum[i + 1]))
  
  step_2_y <- y_coordinate_sum[i + 1] : y_coordinate_sum[i + 2]
  step_2_x <- rep(x_coordinate_sum[i + 1], length(y_coordinate_sum[i + 1] : y_coordinate_sum[i + 2]))
  
  step_x <- c(step_x, c(step_1_x[1:length(step_1_x) - 1], step_2_x[1:length(step_2_x) - 1]))
  step_y <- c(step_y, c(step_1_y[1:length(step_1_y) - 1], step_2_y[1:length(step_2_y) - 1]))

}

step_df <- data.frame(x_coor = step_x, 
                      y_coor = step_y)
step_df$coordinates <- paste(step_df$x_coor, step_df$y_coor, sep = ",")
step_df$second_repeat <- duplicated(step_df$coordinates)
step_df[step_df$second_repeat == T, ][1, ]

ggplot(data = step_df, aes(x = x_coor, y = y_coor, color = second_repeat)) +
  geom_point()
