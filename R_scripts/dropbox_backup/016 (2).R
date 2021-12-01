f.letters <- factor(letters)
f.letters
levels(f.letters) = c("0", "1", "2", "d", "e", "f", "g", "h", "i", 
                      "j", "k", "l", "m", "n", "o", "p", "q", "r", 
                      "s", "t", "u", "v", "w", "x", "y", "z")
f.letters
str(f.letters)

# str prikazuje strukturu varijable f.letters. Pokazuje broj 
# levela u factoru (26), pokazuje koji su to leveli kao niz 
# charactera i takoğer pokazuje niz brojeva. Kako bi bolje vidio
# što taj niz brojeva znaèi napravio sam donji primjer s manje 
# slova:

letters2 <- c("a", "b", "c", "d", "a", "a", "b", "d")
f.letters2 <- factor(letters2)
f.letters2
str(f.letters2)

# Iz ovog primjera zakljuèujem da leveli u faktoru predstavljaju
# kategorije vrijednosti s kojima ulazimo u faktor. U vektoru 
# letters2 imamo 3 puta slovo "a", a kada vektor pretvorimo u faktor
# pod levels se a pojavljuje samo jednom. Niz brojeva nakon popisa
# levela pokazuje kako su pojedine vrijednosti poredane u faktoru: 
# svaki broj predstavlja broj mjesta u levels na kojem se nalazi
# vrijednost za to mjesto. Npr. ako je levels: "a", "b", "c",
# onda struktura 1 1 2 1 3 znaèi da je faktor: a a b a c:

letters3 <- c("a", "a", "b", "a", "c")
f.letters3 <- factor(letters3)
str(f.letters3)
f.letters3


