mydata <- read.table(file = 'hgnc_complete_set.txt', sep = '\t', 
                     header = T, fill = TRUE, quote = "", 
                     stringsAsFactors = FALSE)

mydata <- mydata[, -4:-17]
mydata <- mydata[, -6:-26]
a <- grepl('withdrawn', mydata$Approved.Symbol)
mydata <- subset(mydata, subset = !a)

My.b.s.tree.add <- function(parent, newnode){
  if (parent$Approved.Symbol == newnode$Approved.Symbol){
    return(parent)
  }
  else {
    if (parent$Approved.Symbol > newnode$Approved.Symbol){
      if (length(parent$left) == 0)
      {
        parent$left <- newnode
        return(parent)
      }
      parent$left <- my.b.s.tree.add(parent$left, newnode)
      return(parent)
    }
    else {
      if (length(parent$right) == 0){
        parent$right <- newnode
        return(parent)
      }
      parent$right <- My.b.s.tree.add(parent$right, newnode)
      return(parent)
    }
  }
}



mydata1 <- mydata
mydata1 <- mydata1[order(mydata1$Approved.Symbol), ]

mydata1 <- mydata1[-11:-39308, ]

my.parent <- as.list(my.tree)
my.newnod <- as.list(mydata1[10, ])

my.tree <- My.b.s.tree.add(my.parent, my.newnod)

