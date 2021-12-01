a <- 1:5
b <- 1:3
d <- c(T,F)


a+b
# rezultat je 2 4 6 5 7
# program zbroji vektore, ali uz upozorenje da duljina vektora a nije višekratnik duljine vektora b
# postupak zbrajanja:
#    a  1 2 3 4 5 
#    b  1 2 3 1 2 <- vektor b je kraæi od vektora a pa se prva dva èlana vektora b ponavljaju
#    =  2 4 6 5 7


b+d
# rezultat je 2 2 4
# program zbroji vektore, ali uz upozorenje da duljina vektora b nije višekratnik duljine vektora d
# pri zbrajanju program logièkim podacima iz vektora d pridružuje 1 ako su TRUE ili 0 ako su FALSE
# postupak zbrajanja:
#    b  1 2 3
#    d  T F T <- vektor d je kraæi od vektora b pa se prvi èlan vektora d ponavlja
#       1 0 1
#    =  2 2 4


a+d
# rezultat je 2 2 4 4 6
# program zbroji vektore, ali uz upozorenje da duljina vektora a nije višekratnik duljine vektora d
# pri zbrajanju program logièkim podacima iz vektora d pridružuje 1 ako su TRUE ili 0 ako su FALSE
# postupak zbrajanja:
#    a  1 2 3 4 5
#    d  T F T F T <- vektor d je kraæi od vektora b pa se vektor d ponavlja 2.5 puta
#       1 0 1 0 1
#    =  2 2 4 4 6
  