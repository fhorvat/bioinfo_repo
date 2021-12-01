# Bar-plot se sastoji od stupaca na grafu koji se nalaze iznad 
# oznaka koje predstavljaju kvalitativne varijable. Visina svakog 
# stupca predstavlja velièinu skupine definirane oznakom. Funkcija
# za crtanje je "barplot".

# Histogram se takoðer sastoji od stupaca na grafu ispod kojih
# su oznake koje oznaèavaju kvantitativne varijable. Oznaka ispod 
# stupca može predstavljati jednu vrijednost ili raspon (range) 
# vrijednosti. Visina stupca predstavlja velièinu grupe definiranu
# oznakom ispod stupca. Funkcija za crtanje je "hist". 

# Razlika je što stupac u bar-plotu predstavlja kvalitativnu
# varijablu, a stupac u histogramu predstavlja kvantitativnu 
# varijablu. 

a <- c(1,3,3,2,4,5,5,6,7,10)
barplot(a, main="barplot", names.arg=c("a","b","c","d","e","f","g","h","i","j"))
hist(a, main="histogram")
