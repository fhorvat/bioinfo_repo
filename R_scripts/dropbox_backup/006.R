library(dplyr)
# A)
chr21 <- read.table("002.gtf", stringsAsFactors = F, sep = "\t")

# razdvajanje gene_id, transcript_id stupca na dva odvojena stupca
colnames(chr21) <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
tmp <- matrix(unlist(strsplit(as.character(chr21$"i"), '; ')),
              ncol = 2, byrow = TRUE)
chr21 <- cbind(chr21, as.data.frame(tmp))
chr21 <- select(chr21, -9)
colnames(chr21)[9:10] <- c("gene_id", "transcript_id")

# brisanje "gene_id" i "transcript_id" iz svakog reda u stupcima gene_id i transcript_id
chr21 <- chr21 %>%
  mutate(gene_id = gsub("gene_id ", "", x = chr21$"gene_id")) %>%
  mutate(transcript_id = gsub("transcript_id ", "", x = chr21$"transcript_id"))

# B) 
gene.transc <- read.table("mart_export (1).txt", stringsAsFactors=F, sep="\t",
                           header = TRUE)

# C) All genes that are protein coding
proteins.cod.genes <- gene.transc %>%
  filter(Gene.Biotype == "protein_coding") %>%
  group_by(Ensembl.Gene.ID) %>%
  summarize()





