x <- 1:10
sample(x, size=5)

# Funkcija "sample" uzima nasumièan uzorak velièine n iz vektora 
# duljine m. Ako ne preciziramo "size" uzima nasumièan uzorak duljine
# jednake duljini vektora (permitira vektor): 

x <- 1:10
sample(x)

# Parametar "replace" dozvoljava da jedan èlan vektora bude izabran
# više puta:

x <- 1:10
sample(x, replace=TRUE)

#Pomoæu tog vektora moguæe je uzeti nasumièan uzorak 
# vektora koji je veæi od same velièine vektora:

x <- 1:10
sample(x, size=15, replace=TRUE)
