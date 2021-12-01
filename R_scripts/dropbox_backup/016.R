lst <- list(vct1 = 1:5, vct2 = 5:73, vct3 = seq(from = - pi, to = exp(1), length.out = 24))

mean.lapply <- lapply(lst, mean)
mean.lapply

mean.sapply <- sapply(lst, mean)
mean.sapply

  # lapply vraæa rezultat u obliku liste na kojoj su vrijednosti
  # srednjih vrijednosti za pojedinaène vektore s liste lst

  # sapply je pojednostavljena funkcija koja vraæa vrijednosti
  # srednjih vrijednosti za pojedinaène vektore s liste lst
  # u obliku vektora s imenima vektora s liste lst

mean.lst.for <- list()
for (i in 1:length(lst)){
  mean.lst.for[i] <- mean(lst[[i]])
}
mean.lst.for

x <- 1
mean.lst.while <- list()
while (x < length(lst)+1){
  mean.lst.while[x] <- mean(lst[[x]])
  x <- x + 1
}
mean.lst.while
