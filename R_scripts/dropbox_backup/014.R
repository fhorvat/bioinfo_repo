my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)

write.table(my.data.Char, file = "my.data.Char1.txt",
            sep = ";")

write.table(my.data.Char, file = "my.data.Char2.txt",
            sep = "\t")
