### INFO: get expression of all genes in lncKO data
### DATE: Wed May 23 19:20:31 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/expression")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(tibble)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(Rsamtools)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set input path
inpath <- getwd()

# set output path
outpath <- getwd()

# gtf path
gtf_path <- list.files(path = genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.gtf.gz$", full.names = T)

# sample table path
sample_table_path <- file.path(inpath, "lncKO.all.sample_table.csv")

# summarizedOverlaps path 
se_path <- file.path(inpath, "ensembl.91.GRCm38.p5.20180512.lncKO.summarizedOverlaps.transcripts.RDS")

######################################################## READ DATA
# read gtf
gtf_df <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 

# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

# read summarizedOverlaps
se <- readRDS(se_path)

######################################################## MAIN CODE
### take transcript with most exons for each gene
# get exons for each transcript
gtf_gr <- gtfToGRanges(gtf_df, "exon")

# get unique gene_id/transcript_id combinations
gids <- unique(values(gtf_gr)[c('gene_id', 'transcript_id')])

# splits exon ranges on transcripts
gtf_gr <- split(gtf_gr, gtf_gr$transcript_id)

# orders transcripts based on number of exons in transcripts
gtf_gr <- gtf_gr[order(elementNROWS(gtf_gr), decreasing = T)] 

# keeps only first transcript of each gene (the one with most exons)
gtf_gr <- gtf_gr[!duplicated(gids$gene_id[match(names(gtf_gr), gids$transcript_id)])]

# ### summarizeOverlaps
# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# names(bamfiles) <- sample_table$sample_id
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
# se <- GenomicAlignments::summarizeOverlaps(features = gtf_gr,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = TRUE,
#                                            ignore.strand = TRUE)
# 
# # save RDS
# saveRDS(se, file = file.path(outpath, "ensembl.91.GRCm38.p5.20180512.lncKO.summarizedOverlaps.transcripts.RDS"))

### prepare data
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size) %>% 
  # dplyr::mutate(genotype = ifelse(genotype == "WT", str_c(genotype, "_", str_extract(sample_id, "201[6,7]{1}[:alpha:]{3}")), genotype)) %>%
  dplyr::filter(str_detect(sample_id, "2016Nov")) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# get gene_id of protein coding genes
protein_genes <- 
  gene_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .)

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
se_filt <- se_filt[, match(rownames(sample_table_dds), colnames(se_filt))]
colData(se_filt) <- DataFrame(sample_table_dds)


### DESeq2
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~genotype)

# run DESeq
dds_deseq <- DESeq(dds)

## get results
# set results to fetch
comparison <- c("Lnc21Null", "WT")

# get results, shrink logFC
dds_shrink <- lfcShrink(dds_deseq, contrast = c("genotype", comparison))

# MA plot
# get table for plot
plot_df <-
  dds_shrink %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  tibble::as.tibble(.) %>% 
  dplyr::left_join(., gene_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>% 
  dplyr::arrange(padj) %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_name) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                sign = ifelse(padj < 0.1, "yes", "no"),
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, sign == "no", "not_sign"), 
                regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>% 
  dplyr::arrange(regulation)

# get labels
labels_df <- 
  plot_df %>% 
  dplyr::filter(sign == "yes")

# plot
ggplot() + 
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation), size = 3, shape = 20) +
  geom_label_repel(data = labels_df, aes(x = mean, y = lfc, label = gene_name), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  scale_x_log10(limits = c(1e-01, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-7, 5),
                     breaks = c(-7:5)) +
  scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
  guides(color = FALSE) +
  ggtitle(str_c(comparison[1], " vs. ", comparison[2])) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggsave(filename = file.path(outpath, str_c("MAplot.", comparison[1], "_vs_", comparison[2], ".trans.png")), width = 10, height = 10)


