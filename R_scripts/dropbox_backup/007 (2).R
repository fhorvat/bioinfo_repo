my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)
a <- grep("\\d$", my.data.Char$Approved.Symbol, value = T)
a
  # a se sastoji do svih approved symbols koji na zadnje mjestu imaju
  # znamenku


b <- grep("\\d.$", my.data.Char$Approved.Symbol, value = T)
b
  # b se sastoji do svih approved symbols koji na predzadnjem mjestu 
  # imaju znamenku

  # value je parametar funkcije grep koji ako je TRUE vraæa vrijednost
  # koja odgovara zadanom patternu. Ako je FALSE vraæa
  # mjesto (lokaciju, index) na kojoj se nalazi vrijednost koja 
  # odgovara patternu. Pretpostavljena vrijednost je FALSE. Primjer:

xyz <- c("plavo", "zeleno", "crveno", "uto", "bijelo")
x <- grep("l", xyz, value = TRUE)
y <- grep("l", xyz, value = FALSE)
x
y

  # x sadri sve rijeèi iz character vektora xyz koje sadre slovo l
  # y sadri broj mjesta (smještaj, index) u character vektoru xyz 
  # na kojem se nalazi rijeè koja sadri slovo l 