library("dplyr")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/Fugaku")
rm(list = ls())

# read fpkm data
fpkm_df <- read.csv("fpkm_fugaku_all_samples.csv", stringsAsFactors = F, row.names = 1)
fpkm_df$gene_id <- rownames(fpkm_df)
rownames(fpkm_df) <- NULL
fpkm_df$gene_symbol <- mapIds(org.Mm.eg.db,
                              keys = fpkm_df$gene_id,
                              column ="SYMBOL",
                              keytype = "ENTREZID",
                              multiVals = "first")

# filtering stages
fpkm_df_filtered <- 
  fpkm_df %>%
  dplyr::select(matches(".GV.WE$|.MII.WE$|.1cell.WE$|s_1cell.WE_DNAm"), gene_symbol, gene_id) 
  # mutate(maternal = fpkm_df %>% dplyr::select(matches("GV|MII")) %>% rowMeans(),
  #        zygotic = fpkm_df %>% dplyr::select(matches("2cell|4cell")) %>% rowMeans(),
  #        embrional = fpkm_df %>% dplyr::select(matches("Molura|Blast")) %>% rowMeans())

# fpkm_df_filtered$max_expression <- 
#   fpkm_df_filtered %>% 
#   dplyr::select(1:4) %>% apply(1, which.max)

########################################################################## filtering
# keep genes which are:
# - max. 1C
# - at least 5X higher in 1C than all other

# filtering based on criteria above
fpkm_df_1cell <- 
  fpkm_df_filtered %>%
  filter(s_1cell.WE > (5 * s_1cell.WE_DNAm), 
         s_1cell.WE > (5 * s_GV.WE), 
         s_1cell.WE > (5 * s_MII.WE))

########################################################################## getting intron-less genes
# genes without introns
exonsByKnownGenes <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, "gene")
knownGenes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene, single.strand.genes.only = T)

# get data.frames with gene and exons widths, join, filter only intronless genes
knownGenes_df <- data.frame(gene_id = names(knownGenes), gene_width = width(knownGenes))

exons_df <- data.frame(exons_width = sapply(width(exonsByKnownGenes), sum))
exons_df$gene_id <- rownames(exons_df)
rownames(exons_df) <- NULL
exons_df <- exons_df[, c(2, 1)]

knownGenes_intronless <- 
  left_join(knownGenes_df, exons_df, by = "gene_id") %>%
  filter(gene_width == exons_width)

########################################################################## 
# join genes expressed in 1cell with ranges
genes_1cell <- 
  left_join(fpkm_df_1cell, as.data.frame(knownGenes), by = "gene_id") %>%
  dplyr::select(gene_id, gene_symbol, seqnames, start, end, strand,
                s_GV.WE, s_MII.WE, s_1cell.WE, s_1cell.WE_DNAm) %>%
  mutate(gene_name = mapIds(org.Mm.eg.db,
                            keys = intronless_1cell$gene_id,
                            column ="GENENAME",
                            keytype = "ENTREZID",
                            multiVals = "first")) %>%
  mutate(intronless = ifelse(gene_id %in% knownGenes_intronless$gene_id, T, F))

# write table
write.table(genes_1cell, "intronless_1cell_knownGenes_Fugaku.tsv", quote = F, row.names = F, sep = "\t")
write.table(genes_1cell[genes_1cell$intronless == T, ], "intronless_1cell_knownGenes_Fugaku_2.tsv", quote = F, row.names = F, sep = "\t")

