### INFO: counts reads over mature miRNAs (from miRbase .gff)
### DATE: Thu Aug 16 17:26:28 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/miRNA_expression/sample_correlation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)

library(GGally)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/miRNA_expression"

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRbase path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

######################################################## READ DATA
# read miRbase gtf
mirbase_gr <- rtracklayer::import.gff(con = mirbase_path) 

######################################################## MAIN CODE
# get ranges of mature miRNA
mirna_gr <- mirbase_gr[mcols(mirbase_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name", "Derives_from")]
mcols(mirna_gr)$unique_name <- str_c(mcols(mirna_gr)$Name, ".", mcols(mirna_gr)$Derives_from)

# create data.frame with miRNA coordinates
mirna_df <- 
  mirna_gr %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(unique_name, coordinates, strand) %>% 
  dplyr::filter(!duplicated(unique_name))

### read summarizedOverlaps .RDS for all experiments and writes results table
# ESC bams
esc_bams <- c("s_ES_Dcr_Rsc_tran_r2.SE.mis_0", "s_ES_Dcr_Rsc_utran_r2.SE.mis_0",
              "s_ES_DcrSom_Rsc_tran_r1.SE.mis_0", "s_ES_DcrSom_Rsc_utran_r1.SE.mis_0", 
              "s_ES_WT_tran_r1.SE.mis_0")

esc_bams <- c("s_ES_Dcr_KO_tran_r1.SE.mis_0", "s_ES_Dcr_KO_utran_r1.SE.mis_0", 
              "s_ES_Dcr_Rsc_tran_r1.SE.mis_0", "s_ES_Dcr_Rsc_tran_r2.SE.mis_0", "s_ES_Dcr_Rsc_tran_r3.SE.mis_0", 
              "s_ES_Dcr_Rsc_utran_r1.SE.mis_0", "s_ES_Dcr_Rsc_utran_r2.SE.mis_0", "s_ES_Dcr_Rsc_utran_r3.SE.mis_0",
              "s_ES_DcrSom_Rsc_tran_r1.SE.mis_0", "s_ES_DcrSom_Rsc_tran_r2.SE.mis_0", "s_ES_DcrSom_Rsc_utran_r1.SE.mis_0","s_ES_DcrSom_Rsc_utran_r2.SE.mis_0",
              "s_ES_WT_tran_r1.SE.mis_0", "s_ES_WT_utran_r1.SE.mis_0")
              

# experiment paths
experiment_paths <- 
  c("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){
  
  experiment <- "ES_DcrTrans_2012"
  
  ### read data
  # set summarizedOverlaps .RDS path
  se_path <- file.path(inpath, str_c("miRbase.", experiment, ".se.RDS"))
  
  # library size path
  library_size_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/MosIR_expression", 
                                 str_c("mosIR.", experiment, ".counts_summary.xlsx"))
    
  # sample table path
  sample_table_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD",
                                 experiment, "Data/Documentation", str_c(experiment, ".sample_table.csv"))
  
  # read summarizeOverlaps .RDS
  se <- readRDS(file = se_path)
  
  # read library size df
  library_size_df <- 
    xlsx::read.xlsx(file = library_size_path, sheetName = "library_sizes") %>% 
    tibble::as.tibble(.) %>% 
    tidyr::gather(library_type, library_size, -sample_id) %>% 
    tidyr::separate(library_type, into = c("mis", "library_type"), sep = "\\.") %>% 
    tidyr::unite(sample_id, sample_id, mis, sep = ".") %>% 
    tidyr::spread(library_type, library_size)
  
  # read sample table
  sample_table <-
    readr::read_csv(sample_table_path) %>%
    tidyr::unite(genotype, genotype, transfection, sep = ".") %>%
    dplyr::select(sample_id, genotype)
  
  
  ### calculate RPM
  # get data.frame of counts
  mirna_rpm <- 
    assay(se) %>% 
    as.tibble(.) %>% 
    magrittr::set_colnames(str_remove(colnames(.), ".bam")) %>% 
    dplyr::mutate(mirna_id = mcols(mirna_gr)$unique_name) %>% 
    dplyr::select(mirna_id, everything()) %>% 
    dplyr::filter(!duplicated(mirna_id)) %>% 
    tidyr::gather(sample_id, count, -mirna_id) %>% 
    dplyr::filter(str_detect(sample_id, "mis_0")) %>%
    dplyr::left_join(., library_size_df, by = "sample_id") %>%
    dplyr::mutate(library_size = (`21to23nt_reads` / 1e6), 
                  rpm = count / library_size) %>% 
    dplyr::select(mirna_id, sample_id, rpm) %>% 
    tidyr::spread(sample_id, rpm) %>%
    dplyr::left_join(., mirna_df %>% dplyr::select(-strand), by = c("mirna_id" = "unique_name")) %>% 
    dplyr::select(-c(mirna_id, coordinates)) 
  
  # filter bams in ESC
  if(experiment == "ES_DcrTrans_2012"){
    
    # filter samples
    mirna_rpm %<>%
      dplyr::select(esc_bams)
    
    # change genotype for colors in PCA
    sample_table %<>% 
      dplyr::mutate(genotype = str_remove(genotype, "\\..*$"))
    
  }
  
  ### correlation matrix plot
  # create plot
  cor_pairs <- GGally::ggpairs(mirna_rpm %>% magrittr::set_colnames(., str_remove_all(colnames(.), ".SE.mis_0|^s_T3T_|^s_ES_|r(?=[1-9]{1})") %>% 
                                                                      str_replace_all(., "_", " ") %>% 
                                                                      str_replace(., "(?<=tran) ", ". ")),
                               diag = "blank")
  
  # limit axis on all plots
  for(i in 2:cor_pairs$nrow) {
    for(j in 1:(i - 1)) {
      cor_pairs[i, j] <-
        cor_pairs[i, j] +
        scale_x_continuous(limits = c(0, 500)) +
        scale_y_continuous(limits = c(0, 500))
    }
  }
  
  # add themes
  cor_pairs <- 
    cor_pairs + 
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) 
  
  # save plot
  png(filename = file.path(outpath, str_c("corMatrix.", experiment, ".all.miRBase.RPM.png")), width = 10, height = 10, units = "in", res = 300)
  cor_pairs
  dev.off()
  
  
  ### distance heatmap
  # calculate distance
  dist_df <-
    mirna_rpm %>%
    magrittr::set_colnames(., str_remove_all(colnames(.), ".SE.mis_0|^s_T3T_|^s_ES_|r(?=[1-9]{1})") %>% 
                             str_replace_all(., "_", " ") %>% 
                             str_replace(., "(?<=tran) ", ". ")) %>% 
    t(.) %>%
    dist(.)
  
  # make matrix
  dist_matrix <- as.matrix(dist_df)
  colnames(dist_matrix) <- NULL
  
  # plot
  pheatmap::pheatmap(dist_matrix,
                     clustering_distance_rows = dist_df,
                     clustering_distance_cols = dist_df,
                     col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                     filename = file.path(outpath, str_c("distHeatmap.", experiment, ".all.miRBase.RPM.png")),
                     height = 10,
                     width = 12)

  
  ## PCA - RPM
  # calculate pca
  pca <-
    mirna_rpm %>%
    dplyr::mutate_all(funs(log2(. + 1))) %>%
    t(.) %>%
    stats::prcomp(.)
  
  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
  
  # make tabble for plot
  pca_df <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           sample_id = colnames(mirna_rpm) %>% str_remove(., ".mis_0")) %>%
    dplyr::left_join(sample_table , by = "sample_id") %>%
    dplyr::mutate(genotype = str_replace(genotype, "\\.", " ") %>% 
                    str_replace(., "(?<=tran)", "."), 
                  sample_id = str_remove_all(sample_id, ".SE|^s_T3T_|^s_ES_|r(?=[1-9]{1})") %>% 
                    str_replace_all(., "_", " ") %>% 
                    str_replace(., "(?<=tran) ", ". "))
  
  # plot
  ggplot(data = pca_df, aes(x = PC1, y = PC2, label = sample_id, fill = genotype)) +
    geom_point(size = 7.5, shape = 21) +
    geom_label_repel(fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", show.legend = F) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    # guides(fill = FALSE, shape = FALSE) +
    guides(fill = guide_legend(override.aes = list(shape = 23, size = 5)),
           shape = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = - 0.2),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(outpath, str_c("PCA.", experiment, ".all.miRBase.RPM.png")), width = 11, height = 10)

  # ## PCA - rlog
  # # prepare sample table for DESeq colData
  # sample_table_dds <-
  #   sample_table %>%
  #   as.data.frame(.) %>%
  #   set_rownames(., .$sample_id)
  # 
  # # filter se
  # se_filt <- se[, str_detect(colnames(se), "\\.mis_0")]
  # colnames(se_filt) <- str_remove(colnames(se_filt), "\\.mis_0\\.bam")
  # se_filt <- se_filt[, base::match(rownames(sample_table_dds), colnames(se_filt))]
  # 
  # # add sample table to se
  # colData(se_filt) <- DataFrame(sample_table_dds)
  # 
  # # make DESeqDataSet
  # dds <- DESeqDataSet(se_filt, design = ~genotype)
  # 
  # # rlog transformed counts
  # rlog_df <-
  #   rlog(dds, blind = T) %>%
  #   assay(.)
  
}


