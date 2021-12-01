library("igraph")

####################### PART 1
# Exercise 1
# Create an empty directed graph with 5 nodes. Set color of all nodes to yellow and shape to sphere.
g <- make_empty_graph(n = 5, directed = TRUE)
V(g)$color <- "yellow"
V(g)$shape <- "sphere"


# Exercise 2
# Add the following edges to the graph: 1->2, 1->3, 2->4, 3->4, 4->5.
g <- add.edges(g, c(1,2, 1,3, 2,4, 3,4, 4,5))


# Exercise 3
# Add a vertex to the graph, color it to red and add edges: 3->6, 6->5. Set vertex shape to sphere.
g <- add.vertices(g, 1, color = "red", shape = "sphere")
g <- add.edges(g, c(3,6, 6,5))


# Exercise 4
# Replace edge 1->3 with the edge 3->1.
E(g)
## + 7/7 edges:
## [1] 1->2 1->3 2->4 3->4 4->5 3->6 6->5
g <- delete.edges(g, c(2))
g <- add.edges(g, c(3,1))


# Exercise 5
# Name vertices with letters A-F. List all vertices and edges.
V(g)$name <- LETTERS[1:6]
V(g)
## + 6/6 vertices, named:
## [1] A B C D E F
E(g)
## + 7/7 edges (vertex names):
## [1] A->B B->D C->D D->E C->F F->E C->A


# Exercise 6
# What is the shortest path between vertices A and E.
shortest_paths(g, "A", "E", output = "epath")$epath[1]
## [[1]]
## + 3/7 edges (vertex names):
## [1] A->B B->D D->E


# Exercise 7
# Plot the graph. The size of a vertex should depend on the number of incoming edges.
plot(g, 
     layout = layout_nicely, 
     vertex.size = degree(g, V(g), "in")*15 + 15,
     vertex.label.dist = 0.5, 
     edge.arrow.size = 0.5)


# Exercise 8
# Plot the distribution of vertices degrees.
plot(degree_distribution(g), 
     main = "Degree distribution", 
     xlab = "Degree", 
     ylab = "Frequency")


# Exercise 9
# Create a heatmap from the data.
pal <- colorRampPalette(c("lightblue", "blue"))
a <- as.matrix(get.adjacency(g))
heatmap(a, Rowv = NA, Colv = "Rowv", col = pal(100))

# Exercise 10
# Create graph from subset of vertices with 1 or more incoming edges and plot it. Set the vertice shape to green box. The size of a vertice should depend on the number of outcoming edges.
sg <- induced_subgraph(g, V(g)[degree(g, V(g), "in") >= 1])
plot(sg, 
     layout = layout_nicely, 
     vertex.size = degree(sg, V(sg), "out")*10 + 15,
     vertex.color = "green", 
     vertex.shape = "square", 
     vertex.label.dist = 0.5, 
     edge.arrow.size = 0.5)


####################### PART 2
# A number of employees in a factory was interviewed on a question: 
# “Do you like to work with your co-worker?”. 
# Possible answers are 1 for yes and 0 for no. 
# Each employee gave an answer for each other employee thus creating an adjecancy matrix. 

# Exercise 1
# Load the data and create an unweighted directed graph from the adjecancy matrix. Name the nodes as letters A to Y. Set node color to yellow and shape to sphere. Set the edge’s color to gray and arrow size to 0.2.
setwd("C:/Users/fhorvat/Dropbox/Praksa bioinfo/Projekti/test")
d <- read.csv("sociogram-employees-un.csv", header=FALSE)
g <- graph.adjacency(as.matrix(d), mode = "directed")
V(g)$name <- LETTERS[1:NCOL(d)]
V(g)$color <- "yellow"
V(g)$shape <- "sphere"
E(g)$color <- "gray"
E(g)$arrow.size <- 0.2


# Exercise 2
# Plot the graph.
plot(g)


# Exercise 3
# Calculate network diameter and average closeness.
diameter(g)
mean(closeness(g))


# Exercise 4
# Calculate average network betweenness.
mean(betweenness(g))


# Exercise 5
# Calculate network density and average degree.
graph.density(g)
mean(degree(g, mode="all"))


# Exercise 6
# Calculate network reciprocity and average transitivity.
reciprocity(g)
mean(transitivity(g))


# Exercise 7
# Calculate average eccentricity of the vertices. What is the average distance between two nodes?
mean(eccentricity(g))
mean_distance(g)


# Exercise 8
# Find the hubs and plot graph with node’s size according to their hubs index. Which employee is the biggest hub?
hs <- hub.score(g)$vector
plot(g, 
     layout = layout.fruchterman.reingold, 
     vertex.size = hs*25)
which.max(hs)


# Exercise 9
# Find the authorities and plot graph with node’s size according to their authority index. Which employee is the biggest authority?
as <- authority.score(g)$vector
plot(g, layout=layout_nicely, vertex.size=as*20)
which.max(as)

# Exercise 10
# Show the nodes that make diameter. Plot these nodes larger and in red. Plot edges on this path thicker in red.
diameter.nodes <- get.diameter(g)
diameter.nodes
V(g)$size <- 20
V(g)[diameter.nodes]$color <- "red"
V(g)[diameter.nodes]$size <- V(g)[diameter.nodes]$size+10
E(g)$width <- 1
E(g, path = diameter.nodes)$color <- "red"
E(g, path = diameter.nodes)$width <- 2
plot.igraph(g, layout = layout.fruchterman.reingold)
