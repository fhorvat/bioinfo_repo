library(dplyr)
mydata <- read.table(file = 'Kmeans_test.tsv', sep = '\t', header = T)

  # initializing varijabli
  x <- list()
  y <- 0
  n <- 3
  
  # biranje n nasumiènih  centroida
  cent <- mydata[sample(nrow(mydata), n), ]
  
  for (i in 1:(n + 25)){
    
    # izraèunavanje udaljenosti svake toèke iz data framea od centroida
    mydata1 <- mydata
    for (j in 1:n){
      mydata1[ , ncol(mydata)+j] <- apply(mydata1, 1, function(x) sqrt(sum((x - cent[j, ])^2)))
    }
    
    # pridruživanje centorida svakoj toèki pomoæu which.min (gleda
    # koja vrijednost centroida je najbliža toèki, tj. koja je najmanja
    # vrijednost)
    mydata1$centroid <- apply(mydata1[, (ncol(mydata) + 1):(ncol(mydata)+n)], 1, which.min)
    
    # izraèunavanje novih centroida pomoæu funkcija iz paketa dplyr. 
    # Prvo grupiram po postojeæim centroidima, zatim izraèunavam srednju 
    # vrijednost. Dobijem data frame sa srednjim vrijednostima svih kolona 
    # pa brišem kolone koje mi ne trebaju. 
    cent.new <- mydata1 %>%
      group_by(centroid) %>%
      summarise_each(funs(mean))
    cent.new <- cent.new[, -(ncol(mydata)+2):-(ncol(mydata)+2+n)]
    cent.new <- cent.new[-1]
    
    # provjeravanje jesu li novi centroidi jednaki postojeæima
    # (ako jesu prekida se petlja, ako nisu funkcija kreæe ispoèetka
    # s novim vrijednostima centroida)
    x[i] <- list(all.equal(cent, cent.new))
    if (x[i] == "TRUE") break
    cent <- cent.new
  }
  
  # dodavanje kolone s dodijeljenim brojem centroida u data.frame 
  # s vrijednostima centroida i preslagivanje kolona tako da kolona
  # s dodijeljenim brojem centroida bude prva kolona
  cent$centroid <- row.names(cent)
  cent <- cent[, c(ncol(cent), 1:(ncol(cent) - 1))]
  
  # konaèni data set s pridruženim vrijednostima centroida
  mydata1 <- subset(mydata1, select = c(1:ncol(mydata), ncol(mydata1)))
  mydata1 <- cbind(mydata1, cent[match(mydata1$centroid, cent$centroid), 2:ncol(cent)])
  
  # plot - samo ako data ima najviše 2 dimenzije
  if(ncol(mydata) < 3){
    plot(mydata1[1:2], col = mydata1$centroid, pch = 20)
    points(cent[2:3], col = 4, pch = 8, cex = 2)
    
    # Voronoi regije
    for (i in 1:n){
      vor <- subset(mydata1, select = c("x", "y"), centroid == 1)
      vor.hull <- lapply(vor, "[", chull(vor))
      polygon(vor.hull, lty = 2, border = i)
    }
  }


vor <- subset(cent, select = c("x", "y"), centroid == 1 | centroid == 2)
vor.hull <- lapply(vor, "[", chull(vor))
lines(vor.hull)

