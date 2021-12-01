xyz <- c("plavo", "zeleno", "crveno", "žuto", "bijelo")

  # Kao što sam objasnio u prošlom zadatku, grep vraæa vrijednost
  # iz character vektora koja odgovara zadnom patternu ili vraæa 
  # smještaj te vrijednosti u vektoru (ovisno o argumentu value)
  
  # grepl u character vektoru provjerava jednu po jednu vrijednost
  # za zadani pattern i vraæa TRUE ili FALSE za svaku vrijednost 
  # vektora s obzirom na to sadrži li zadani pattern
  # (TRUE ako sadrži, FALSE ako ne sadrži). Primjer:

x <- grep("v", xyz, value = TRUE)
y <- grepl("v", xyz)
x
y

  # U primjeru provjeravamo koja rijeè iz vektora xyz sadrži slovo "v". 
  # x je dobiven grep funkcijom i sadrži sve rijeèi koje u sebi
  # imaju slovo "v". 
  # y je dobiven grepl funkcijom i sadrži logièke vrijednosti po redu
  # za svaku rijeè iz vektora xyz s obzirom na to sadrži li slovo "v" ili ne.
  # Npr. na prvom mjestu u xyz je vrijednost "plavo" koja sadrži "v" 
  # pa je prva vrijednost u logièkom vektoru y TRUE. Na drugom mjestu
  # u xyz je "zeleno" što ne sadrži "l" pa je zato na drugom mjestu u 
  # logièkom vektoru y FALSE itd. za sve vrijednosti iz vektora xyz.
  

