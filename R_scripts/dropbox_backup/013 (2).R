set.seed (2)
a <- rnorm(5)
a

  # set.seed je funkcija koja postavlja "seed" za generator 
  # nasumiènih brojeva (RNG) u R-u. RNG u R-u je ustvari pseudo-RGN
  # koji generira brojeve algoritmom koji za isti "seed" daje 
  # uvijek istu sekvencu brojeva. Tako su nasumièno generirani 
  # brojevi pod istim "seedom" jednaki i ostaju jednaki bez obzira 
  # na duljinu sekvence. 
  
  # Korištenje set.seed funkcije je korisno kada elimo nasimiènu 
  # sekvencu brojeva koja æe biti ponovljiva svaki put.  

set.seed (2)
b <- rnorm(5)
a
b

  # a i b sadre jednake brojeve jer su generirani s istim "seedom"

set.seed (2)
d <- rnorm(7)
d
a

  # d sadri istu sekvencu prvih 5 brojeva kao a i b, a sekvenca se
  # nastavlja i dalje jer je d dui od a i b

set.seed(1)
e <- rnorm(5)
e
a

  # e sadri razlièite brojeve od a, b i d jer je generiran s
  # razlièitim seedom 

