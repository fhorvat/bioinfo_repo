my.data.Char <- read.table("gene_with_protein_product.txt", header = T, 
                           sep =  "\t", fill = TRUE, quote = "", 
                           stringsAsFactors = F)

x <- ifelse(my.data.Char$Previous.Symbols == "", 
            my.data.Char$Approved.Symbol, 
            paste(my.data.Char$Approved.Symbol, ",", 
                  my.data.Char$Previous.Symbols))
y <- grep(" ", x)
x <- substr(x, y, y+1)

my.data.Char$all.symbols <- x



