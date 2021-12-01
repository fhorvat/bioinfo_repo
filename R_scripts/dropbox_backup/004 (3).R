my.data.read.table <- read.table("gene_with_protein_product.txt", header = T,
                                 sep =  "\t", fill = TRUE, quote = "")

my.data.read.csv <- read.csv("gene_with_protein_product.txt", header = T, 
                             sep =  "\t")

my.data.read.delim2 <- read.delim2("gene_with_protein_product.txt")

  # read.csv i read.delim kao oznaku za decimalni znak èitaju
  # ".", a read.csv2 i read.delim2 ","

  # stringsAsFactors je logical argument koji kae R-u treba li
  # stringove iz filea kojeg uèitavamo pretvoriti u faktore
  # (ako je TRUE) ili u charactere (ako je FALSE).
  # Pretpostavljena vrijednost je TRUE. 

my.data.Fact <- read.table("gene_with_protein_product.txt", header = T,
                                 sep =  "\t", fill = TRUE, quote = "")
a <- sapply(my.data.Fact, class)
a

my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)
b <- sapply(my.data.Char, class)
b

