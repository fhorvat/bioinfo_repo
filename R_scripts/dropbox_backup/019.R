Monty.hall <- function(){
  nagrade <- c("koza", "koza", "auto")
  vrata <- sample(nagrade, 3)
  vrata <- vrata[-1]  # jer natjecatelj uvijek bira vrata 1
  x <- match("koza", vrata)  # iza kojih vrata se nalazi koza 
  vrata <- vrata[-x]  # voditelj otvara vrata iza kojih se nalazi koza
  # ostaju samo jedna vrata
  # to su vrata koja natjecatelj može izabrati kada mu je ponuðena zamjena
  return(vrata) 
  
}
  # 100000 puta pokrenemo funkciju Monty.Hall, 
  # rezultati se spremaju:

rezultat <- replicate(100000, Monty.hall()) 
  
  # raèunam kolika je šansa da je natjecatelj pobjedio za oba sluèaja:

  # natjecatelj se odluèuje na zamjenu i pobjeðuje 
  # (što znaèi da je iza "zadnjih preostalih" vrata iz 
  # funkcije bio auto):

zamjena.pobjeda <- rezultat[rezultat == "auto"]
length(zamjena.pobjeda) # broj pobjeda ako se natjecatelj odluèi na zamjenu
vjerojatnost.zamjena.pobjeda <- length(zamjena.pobjeda)/100000
vjerojatnost.zamjena.pobjeda


  # natjecatelj ostaje pri prvom izboru i pobjeðuje
  # (iza prvih vrata je auto, iza "zadnjih" vrata je koza):
prvi.izbor.pobjeda <- rezultat[rezultat == "koza"]
length(prvi.izbor.pobjeda) # broj pobjeda ako natjecatelj ne mijenja vrata
vjerojatnost.prvi.izbor.pobjeda <- length(prvi.izbor.pobjeda)/100000
vjerojatnost.prvi.izbor.pobjeda
