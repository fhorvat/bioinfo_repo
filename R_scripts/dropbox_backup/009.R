my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)

my.data.Char$my.hgnc <- substr(my.data.Char$HGNC.ID, 6, 
                               nchar(my.data.Char$HGNC.ID))

