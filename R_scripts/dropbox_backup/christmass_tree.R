library(ggplot2)

tree_count <- 5000

r1 <- runif(tree_count, 0, 1)
r2 <- runif(tree_count, 0, 1)
x = (sqrt(r1) * (1 - r2)) * 1 + (sqrt(r1) * r2) * 0.5
y = (sqrt(r1) * r2) * 1
data <- data.frame(x = x, y = y, 
                   color = rep("tree",tree_count), 
                   size = rep(1, tree_count))

decoration_colors <- stringr::str_c("decoration", 1:6)
decoration_count <- 175
r1 <- runif(decoration_count, 0, 1)
r2 <- runif(decoration_count, 0, 1)
x = (sqrt(r1) * (1 - r2)) * 1 + (sqrt(r1) * r2) * 0.5
y = (sqrt(r1) * r2) * 1
decoration <- data.frame(x = x, y = y, 
                         color = sample(decoration_colors, decoration_count, replace = T), 
                         size = rep(1, decoration_count))

data <- rbind(data, decoration)
table(data$color)

ggplot(data, aes(x, y, color = color)) +
  geom_point(size = 5) + 
  scale_color_manual(values = c("#3C8D0D","#BE0000", "#BE0000", 
                                "#FA6900", "#F8CA00", "#F8CA00", 
                                "#107FC9"))