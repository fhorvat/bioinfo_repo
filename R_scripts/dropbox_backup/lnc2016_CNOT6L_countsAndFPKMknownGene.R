library(dplyr)
library(readr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
# show_col()

options(bitmapType = 'cairo')
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Analysis")

# files lnc Nov2016 experiment
sample_table_lnc2016 <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/RNAseq_2016_11_23_sampleTable.csv", header = T) %>%
  dplyr::select(ID = sample, Time.Course = Sample, Treatment.Control = Type) %>%
  mutate(Treatment.Control = gsub(" B6", "", Treatment.Control), 
         ID = gsub("_16.*|_[A,C,T,G].*", "", ID), 
         name = paste(gsub(".*_[1-8]_", "", ID), Treatment.Control, sep = "_"), 
         experiment = "lnc2016") %>% 
  left_join(data.frame(sample_paths = list.files("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                                 pattern = "*.bam$",
                                                 recursive = T,
                                                 full.names = T), 
                       logs_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                               pattern = "*Log.final.out", 
                                               recursive = T, 
                                               full.names = T), 
                       stringsAsFactors = F) %>%
              mutate(ID = gsub("^/.*/|_16.*|_[A,C,T,G].*", "", sample_paths)) %>%
              rowwise() %>%
              mutate(lib_size = (as.integer(read.delim(logs_paths, header = F, stringsAsFactors = F)[8, 2]) / 10^6)), 
            by = "ID") %>% 
  mutate(name = ifelse(grepl("HV", ID), paste0(name, "_new"), name))

# files CNOT6L experiment
sample_table_CNOT6L <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_sample_list_11919R_2015_10_29.csv", header = T) %>%
  dplyr::select(ID, Time.Course, Treatment.Control) %>%
  mutate(name = paste(ID, Time.Course, Treatment.Control, sep = "_"),  
         experiment = "CNOT6L") %>%
  left_join(data.frame(sample_paths = list.files("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                                 pattern = "*.bam$",
                                                 recursive = T, 
                                                 full.names = T), 
                       logs_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                               pattern = "*Log.final.out", 
                                               recursive = T, 
                                               full.names = T), 
                       stringsAsFactors = F) %>%
              mutate(ID = gsub("^/.*/|_.*", "", sample_paths)) %>%
              rowwise() %>%
              mutate(lib_size = (as.integer(read.delim(logs_paths, header = F, stringsAsFactors = F)[8, 2]) / 10^6)), 
            by = "ID")

# combine_tables
sample_table_all <- rbind(sample_table_lnc2016, sample_table_CNOT6L)

# exons by genes
ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")

################################################################## counting
# counts lnc2016 samples
bamfiles_lnc2016 <- BamFileList(sample_table_all[sample_table_all$experiment == "lnc2016", "sample_paths"], yieldSize = 2000000)
register(MulticoreParam())
se_lnc2016 <- summarizeOverlaps(features = ebg, 
                                reads = bamfiles_lnc2016, 
                                mode = "Union", 
                                singleEnd = TRUE, 
                                ignore.strand = TRUE)

# counts CNOT6L
bamfiles_CNOT6L <- BamFileList(sample_table_all[sample_table_all$experiment == "CNOT6L", "sample_paths"], yieldSize = 2000000)
register(MulticoreParam())
se_CNOT6L <- summarizeOverlaps(features = ebg, 
                               reads = bamfiles_CNOT6L, 
                               mode = "Union", 
                               singleEnd = FALSE, 
                               ignore.strand = TRUE)

################################################################## counts and fpkm tables
# counts data.frame
count_df <- cbind(as.data.frame(assay(se_lnc2016)), as.data.frame(assay(se_CNOT6L)))
colnames(count_df) <- sample_table_all$"ID"

# fpkm data.frame
fpkm_df <- count_df
fpkm_df$width <- sapply(width(ebg), sum)
invisible(lapply(X = as.character(sample_table_all$"ID"), 
                 FUN = function(X) fpkm_df[, X] <<- fpkm_df[, X] / (sample_table_all[sample_table_all$ID == X, "lib_size"] * (fpkm_df$width / 1000))))
colnames(fpkm_df) <- c(sample_table_all$"name", "width")
fpkm_df$entrezID <- rownames(fpkm_df)
rownames(fpkm_df) <- NULL
write.csv(fpkm_df, "fpkm_knownGene_mm10_lnc2016_CNOT6L.csv", quote = F, row.names = F)

colnames(count_df) <- sample_table_all$"name"
count_df$entrezID <- rownames(count_df)
rownames(count_df) <- NULL
write.csv(count_df, "counts_knownGene_mm10_lnc2016_CNOT6L.csv", quote = F, row.names = F)
