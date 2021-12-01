rivers
sum.for <- 0
for (i in 1:length(rivers)){
  if (rivers[i] > 650) {
    sum.for <- sum.for + rivers[i]
  }
  else{
    sum.for
  }
}
sum.for

sum.subset <- sum(subset(rivers, rivers > 650))
sum.subset
