my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)

my.data.Char$my.approved.symbol <- sub("CYP", "XXX", 
                                       my.data.Char$Approved.Symbol)

  # provjera:
x <- my.data.Char[grep("XXX", my.data.Char$my.approved.symbol), ]
x

  # sub zamjenjuje pattern u stringu zadanom zamjenom (replacement)
  # prvi put kad se pattern pojavi u stringu
  
  # gsub zamjenjuje pattern u stringu zadanom zamjenom svaki puta
  # kad se pattern pojavi u stringu. Primjer:

xyz <- c("Ivan", "Ivana", "Anastazija", "Valentina", "Tea")
y <- sub("a", "@", xyz)
z <- gsub("a", "@", xyz)
y
z

  # U y je prvo slovo "a" u imenima zamijeneno s @
  # U z je svako slovo "a" u imenima zamijenjeno s @
 