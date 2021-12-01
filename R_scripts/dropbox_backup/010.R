my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)
a <- my.data.Char[c(grep("RNASELI", my.data.Char$Approved.Symbol),
                    grep("RNASELI", my.data.Char$Previous.Symbols)), ]
a
