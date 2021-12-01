my.hclust <- function(hclust.df){
  
  # Sreðivanje imena redova data.framea, inicijalizacija varijabli
  rownames(hclust.df) <- gsub(" ", "", rownames(hclust.df))
  hclust.df1 <- hclust.df
  hclust.df1.nrow <- nrow(hclust.df1)
  merge.mat <- matrix(c(0, 0), ncol = 2)
  height.vector <- NULL
  order.vector <- NULL
  check.value <- NULL
  dist.mat.vector <- NULL
  
  # Funkcija za raèunanje udaljenosti za prvi distance matrix
  eucl <- function (x, y){
    sqrt(sum((hclust.df1[x, ] - hclust.df1[y, ])^2))
  }
  
  # Prvi distance matrix
  dist.matrix <- as.data.frame(outer(1 : nrow(hclust.df1), 
                                     1 : nrow(hclust.df1), Vectorize(eucl)))
  names(dist.matrix) <- row.names(hclust.df1[1:nrow(hclust.df1), ])                        
  row.names(dist.matrix ) <- names(dist.matrix)
  dist.matrix[dist.matrix == 0] <- NA
  
  # Funkcija za raèunanje udaljenosti za redove u hclust data.frame-u
  eucl1 <- function (x, y){
    sqrt(sum((x - y)^2))
  }
  
  # Stvaranje klustera i spremanje podataka (height, merge) u petlji
  repeat{
    
    # Traženje najmanje udaljenosti
    x <- rowSums(dist.matrix == min(dist.matrix, na.rm = T), na.rm = T) > 0
    min.dist.names <- rownames(dist.matrix[x, ])
    
    # Stvaranje podataka za merge i stavljanje u matrix
    check.names <- which(rownames(dist.matrix) %in% min.dist.names)
    for (i in 1:2){
      if (check.names[i] <= nrow(hclust.df)){
        check.value[i] <- check.names[i]
        check.value[i] <- check.value[i] * (-1)
      }
      if(check.names[i] > nrow(hclust.df)){
        check.value[i] <- check.names[i] - hclust.df1.nrow
      }
    }
    
    if(merge.mat[1, 1] == 0){
      merge.mat <- matrix(check.value, ncol = 2)  
    } else{
      merge.mat <- rbind(merge.mat, check.value)
    }
           
   
    # Stavljanje height podataka u vektor
    height.value <- dist.matrix[rownames(dist.matrix) %in% min.dist.names[1],
                                colnames(dist.matrix) %in% min.dist.names[2]]
    height.vector <- c(height.vector, height.value)
    
    # Novi cluster
    min.dist.data <- subset(hclust.df1, rownames(hclust.df1) %in% min.dist.names)
    new.cluster <- colMeans(min.dist.data)
    
    # Vezanje podataka o novom klusteru sa starim podacima
    hclust.df1 <- rbind(hclust.df1, new.cluster)
    rownames(hclust.df1)[nrow(hclust.df1)] <- paste(rownames(min.dist.data), 
                                                      sep = "", collapse = " ")
    
    # Brisanje starih podataka iz arrest data.framea i distance matrixa
    hclust.df1[rownames(hclust.df1) %in% min.dist.names, ] <- NA
    dist.matrix[rownames(dist.matrix) %in% min.dist.names, ] <- NA
    dist.matrix[, colnames(dist.matrix) %in% min.dist.names] <- NA
  
    # Raèunanje novih udaljenosti
    for (j in 1 : (nrow(hclust.df1) - 1)){
      dist.mat.vector[j] <- eucl1(hclust.df1[j, ], hclust.df1[nrow(hclust.df1), ])  
    }
    
    # Stavljanje novih udaljenosti u matrix
    dist.matrix <- cbind(dist.matrix, dist.mat.vector)
    dist.mat.vector <- c(dist.mat.vector, NA)
    dist.matrix <- rbind(dist.matrix, dist.mat.vector)
    
    names(dist.matrix) <- row.names(hclust.df1[1:nrow(hclust.df1), ])
    row.names(dist.matrix) <- names(dist.matrix)
    dist.matrix[dist.matrix == 0] <- NA
    
    # Prekidanje petlje kada je merge matrica dugaèka n - 1
    if(nrow(merge.mat) == (hclust.df1.nrow - 1)){
      break
    }    
  }
  
  # Stavljanje order podataka u vektor
  order.names <- rownames(dist.matrix)[nrow(dist.matrix)]
  order.names <- unlist(strsplit(order.names, split = " "))
  order.value <- match(order.names, row.names(hclust.df))
  order.value <- order.value[!is.na(order.value)]
  
  # Sreðivanje merge matrice
  rownames(merge.mat) <- NULL
  names(merge.mat) <- NULL
  attr(merge.mat, "dimnames") <- NULL
  
  # Stvaranje liste i plot
  dendr.list <- list()
  dendr.list$merge <- merge.mat
  dendr.list$height <- height.vector
  dendr.list$order <- order.value
  dendr.list$labels <- rownames(hclust.df)
  class(dendr.list) <- "hclust"  
  dendr.list <- as.dendrogram(dendr.list)
  return(plot(dendr.list))  
}

my.hclust(USArrests)

