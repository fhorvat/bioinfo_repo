library(ggplot2)
options(bitmapType = 'cairo')

rptmsk_viz <- read.delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz", stringsAsFactors = F, header = T)
rptmsk_viz <- rptmsk_viz[, c("genoName", "genoStart", "genoEnd", "strand", "repName", "id")]
colnames(rptmsk_viz) <- c("seqnames", "start", "end", "strand", "element_name", "VIZ_ID")

# filter for MT, ORR and MLT
rptmsk_viz_filtered <- rptmsk_viz[grep("MT|ORR|MLT", rptmsk_viz$element_name), ]
rptmsk_viz_filtered <- rptmsk_viz_filtered[-grep("MMTV-int|MLTR", rptmsk_viz_filtered$element_name), ]

################################################################## split int from LTRs, check if whole insert with -int also contains at least one LTR 
# split int from LTRs
rptmsk_viz_filtered_int <- rptmsk_viz_filtered[grep("-int", rptmsk_viz_filtered$element_name), ]
rptmsk_viz_filtered_ltr <- rptmsk_viz_filtered[-grep("-int", rptmsk_viz_filtered$element_name), ]

# import insert positions
inserts <- rptmsk_viz_filtered_ltr
# inserts$insert_class <- substr(inserts$element_name_paste, 0, 3)
# inserts$insert_class <- gsub("MT[A|B|C|D|E]", "MT", inserts$insert_class)
inserts$insert_class <- "MaLR"
inserts <- inserts[, c("seqnames", "start", "end", "insert_class")]

# import genes position
genes <- as.data.frame(knownGenes_genes)
genes$insert_class <- "genes"
genes <- genes[, c("seqnames", "start", "end", "insert_class")]

# combine into one df
inserts_genes <- rbind(inserts, genes)

# filter by chomosome
chr <- "chr12"
inserts_genes_filtered <- inserts_genes[inserts_genes$seqnames == chr, ]

# make a density plot of genes over the provided chromosomes (or scaffolds ...)
plottedGenes <- 
  ggplot() + 
  geom_histogram(data = inserts_genes_filtered, aes(x = start, fill = insert_class), binwidth = 1000000) + 
  scale_fill_manual(values = c("red", "orange")) +
  ggtitle("Insert density over mm10 chr12") + 
  xlab("Genomic position (bins 1 Mb)") + 
  ylab("Number of inserts")

# save it to an image
png("chr12_insert_genes_density_repeatMaskerVIZ_ltr.png", width = 1500, height = 300)
print(plottedGenes)
dev.off()

