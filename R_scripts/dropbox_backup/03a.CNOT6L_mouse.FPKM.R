### INFO: get expression of all genes in CNOT6L data
### DATE: Sat Mar 10 19:01:11 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# reduced exons path
exons_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.reducedExons.RDS"

# stats and tracks, bam paths
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10_noMultimapFilter"

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read stats and tracks table
stats_df <- readr::read_csv(file = list.files(mapped_path, pattern = "*stats_and_tracks.csv", full.names = T))

# read info about chosen genes
genes_info <- readr::read_csv(file.path(outpath, "CNOT.GRCm38.89.20180305.geneInfo.csv"))

######################################################## MAIN CODE
# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .)

# construct sample table
sample_table <-
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.bam$", full.names = T)) %>%
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".Aligned.sortedByCoord.out.bam", "")) %>%
  dplyr::left_join(stats_df %>% dplyr::select(sample_id, library_size = mapped_genome_except_rDNA), by = "sample_id") %>%
  dplyr::mutate(stage = str_extract(sample_id, "1C|MII|GV"),
                genotype = str_extract(sample_id, "WT|KO")) %>%
  dplyr::select(sample_id, stage, genotype, library_size, bam_path) %T>%
  readr::write_csv(., path = file.path(outpath, "CNOT6L.sample_table.csv"))

# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = FALSE,
#                                            ignore.strand = TRUE)
# saveRDS(se, file = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.RDS"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.RDS"))

# read sample table
sample_table <- readr::read_csv(file = file.path(outpath, "CNOT6L.sample_table.csv"))

### FPKM
# get data.frame of counts, transform to FPKM, write
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".Aligned.sortedByCoord.out.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  dplyr::mutate(stage = str_replace_all(sample_id, "s_|_r[1-3]{1}", "")) %>%
  dplyr::group_by(gene_id, stage) %>%
  dplyr::summarise(average_fpkm = mean(fpkm)) %>%
  dplyr::ungroup(.) %>%
  tidyr::spread(key = stage, value = average_fpkm) %T>%
  readr::write_csv(., path = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"))

# get FPKM of chosen genes, write
fpkm_df %>%
  dplyr::right_join(., genes_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%
  dplyr::select(gene_name, everything(), -gene_id) %>%
  dplyr::arrange(gene_name) %T>%
  readr::write_csv(., path = file.path(outpath, "CNOT.GRCm38.89.CNOT6L.avgFPKM.csv"))

### correlation between replicates
# data for corelation
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".Aligned.sortedByCoord.out.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  tidyr::spread(key = sample_id, value = fpkm) %>% 
  dplyr::select(-gene_id) %>% 
  magrittr::set_colnames(., str_replace(colnames(.), "s_|r", "") %>% str_replace_all(., "_", " "))

# plot correlation matrix for each stage
for(filt_stage in c("1C", "MII", "GV")){
  
  # filter FPKM, limit
  fpkm_df_filt <- 
    fpkm_df %>% 
    dplyr::select(which(str_detect(colnames(.), filt_stage)))
  
  # create plot
  cor_pairs <- GGally::ggpairs(fpkm_df_filt, diag = "blank")
  
  # modify plots
  cor_pairs2 <- cor_pairs
  for(i in 2:cor_pairs2$nrow) {
    for(j in 1:(i-1)) {
      cor_pairs2[i,j] <- 
        cor_pairs2[i,j] +
        scale_x_continuous(limits = c(0, 200)) +
        scale_y_continuous(limits = c(0, 200))
    }
  }
  
  # add themes
  cor_pairs2 <- 
    cor_pairs2 + 
    theme_bw() +
    theme(
      # axis.text.x = element_blank(), 
      # axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) 
  
  # save plot
  ggsave(filename = file.path(outpath, "results", str_c("corMatrix.CNOT6L.", filt_stage, ".GRCm38.89.FPKM.png")), 
         plot = cor_pairs2, width = 10, height = 10)
  
}
