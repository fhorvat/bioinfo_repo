  # Faktor je varijabla u koju su vrijednosti spremljene kao 
  # niz integera sa pridruženim setom character vrijednosti
  # koje se koriste kad se faktor prikazuje (levels). 
  # U faktoru je svaka jedinstvena vrijednost spremljena samo
  # jednom (kategorija), a niz integera pokazuje kako su jedinstvene
  # vrijednosti poredane u faktoru. 

data = c("a", "b", "c", "d", "e", "f", "a", "b", "g", "e", "a")
f.data = factor(data)
f.data

data2 = c("a", "b", "c", "d", "g", "e", "a")
f.data2 = factor(data2)
f.data2

  # U vektor su vrijednosti spremljene po redu kako ih unesemo, 
  # nema kategorija nego je svaka vrijednost spremljena baš onako
  # kako je unesemo u vektor. Vektor je niz podataka indeksiranih
  # po poziciji u jednoj dimenziji. Osim toga, na vektorima možemo
  # koristiti aritmetièke funkcije, a na faktorima ne (jer su
  # faktori kategorizirani). 