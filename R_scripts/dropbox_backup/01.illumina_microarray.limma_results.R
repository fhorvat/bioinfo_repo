### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/Freimer_microarray")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GEOquery)
library(biomaRt)
library(beadarray)
library(beadarrayExampleData)
library(Biobase)
library(limma)
library(Biostrings)
library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### experiment
# set experiment
experiment <- "Freimer_2017_CurrBiol_GSE92658"

# base experiment path
base_experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays"

# experiment path 
experiment_path <- file.path(base_experiment_path, experiment)

# list normalized .csv files
norm_files_list <- list.files(path = experiment_path, pattern = "norm\\.csv\\.gz", full.names = T)

# chip annotation file path
chip_annotation_path <- file.path(experiment_path, "GPL6885.annot.gz")


## gene annotation
# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# transcripts info path
transcripts_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.transcriptInfo.csv$"), full.names = T)

######################################################## READ DATA
# read .csv files
norm_tb <- purrr::map(norm_files_list, function(path){
  
  # read file
  tb <- 
    readr::read_csv(path) %>% 
    dplyr::select(-contains("Pval")) %>% 
    dplyr::rename(probe = X1) %>% 
    as.data.frame(.) %>% 
    tibble::column_to_rownames(., var = "probe") %>% 
    as.matrix(.)
  
}) %>% 
  set_names(c("ctrl", "miR-15a"))

# read chip annotation
chip_annotation <- readr::read_delim(chip_annotation_path, delim = "\t", skip = 28)

# read genes info
genes_info <-
  readr::read_csv(genes_info_path) %>%
  tidyr::unite(., coordinates, seqnames, start, end, sep = " ")

# read transcripts info
transcripts_info <- readr::read_csv(transcripts_info_path)

######################################################## MAIN CODE
### tidy files
# get connection between probe ID's and gene symbols
chip_tidy <- 
  chip_annotation %>% 
  dplyr::select(ID, gene_name = `Gene symbol`, gene_description = `Gene title`, coordinates = `Chromosome annotation`) %>% 
  tidyr::separate(coordinates, into = c("coordinates", "coordinates2"), sep = "///") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("coordinates")), 
                   .funs = list(~str_replace(., "Chromosome ", "chr") %>% 
                                  str_replace(",.+ \\(", " ") %>% 
                                  str_replace("\\.\\.", " ") %>% 
                                  str_remove_all("\\)|, complement"))) %>% 
  dplyr::mutate(gene_name = replace(gene_name, gene_name == "Gucy1a3", "Gucy1a1"))

### experiment
# GEO accession
geo_accession <- str_remove(experiment, ".+_")

# download data
gset <- getGEO(geo_accession, 
               GSEMatrix = TRUE, 
               getGPL = T, 
               AnnotGPL = T, 
               destdir = outpath)

# tidy sample table
sample_table <- 
  gset %>% 
  .[[1]] %>% 
  pData(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(sample_id = title, geo_accession, 
                genotype = characteristics_ch1, 
                treatment = characteristics_ch1.1) %>% 
  dplyr::mutate(stage = "GV", 
                genotype = stringr::str_to_title(genotype) %>% str_replace(., ": ", "_"), 
                treatment = str_remove(treatment, "mir: ") %>% 
                  replace(., . == "none", "ctrl")) %>% 
  dplyr::select(sample_id, geo_accession, stage, genotype, treatment)


### get differentially expressed genes
results_list <- purrr::map(names(norm_tb), function(result){
  
  # filter sample table
  sample_table_filt <-
    sample_table %>%
    dplyr::filter(treatment == result)
  
  # get normalized expression, change column order to match sample table
  exprs_tb <- 
    norm_tb[[result]] %>% 
    .[, match(sample_table_filt$sample_id, colnames(.))]
  
  # check
  if(all(colnames(exprs_tb) != sample_table_filt$sample_id)){
    stop("Order of columns doesn't match order of samples in sample table")
  }
  
  # create design matrix
  f <- factor(sample_table_filt$genotype, levels = c("Cre_Pos", "Cre_Neg"))
  design <- model.matrix(~0 + f)
  colnames(design) <- c("Cre_Pos", "Cre_Neg")
  
  # fit the model
  fit <- lmFit(exprs_tb, design)
  contrast.matrix <- makeContrasts(Cre_PosvsCre_Neg = Cre_Pos - Cre_Neg, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, trend = T, robust = T)
  
  # get the results, tidy
  results_tb <- 
    topTable(fit2, n = nrow(exprs_tb)) %>% 
    as_tibble(., rownames = "ID") %>% 
    dplyr::select(ID, log2FoldChange = logFC, log2_mean_exp = AveExpr, padj = adj.P.Val) %>% 
    dplyr::arrange(padj) %>% 
    dplyr::left_join(., chip_tidy %>% dplyr::select(ID, gene_name), by = "ID") %>% 
    dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_name, gene_description, coordinates), by = "gene_name") %T>%
    readr::write_csv(., file.path(outpath, 
                                  str_c("Freimer_2017", "diff_exp", "CrePos_vs_CreNeg", result, "csv", sep = ".")))
  
  # return
  return(results_tb)
  
}) %>% 
  set_names(names(norm_tb))


### get significant results
results_significant_list <- invisible(purrr::map(names(results_list), function(result){
  
  # filter table
  results_df_sign <-
    results_list[[result]]  %>%
    dplyr::filter(padj <= 0.05) %T>%
    readr::write_csv(., file.path(outpath, 
                                  str_c("Freimer_2017", "diff_exp", "significant", "CrePos_vs_CreNeg", result, "csv", sep = ".")))
  
  # return
  return(results_df_sign)
  
})) %>%
  set_names(., names(results_list))


### plot MA plots
# loop through results
invisible(purrr::map(names(results_list), function(result){
  
  result <- "miR-15a"
  
  ## prepare results
  # get results table
  results_df <- results_list[[result]]
  
  # significant results
  results_df_sign <-
    results_significant_list[[result]] %>%
    dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down")) %>%
    dplyr::select(ID, regulation)
  
  
  ## MA plot
  # data for plot
  plot_df <-
    results_df %>%
    dplyr::select(mean = log2_mean_exp, lfc = log2FoldChange, padj, ID, gene_name) %>%
    dplyr::left_join(., results_df_sign, by = "ID") %>%
    dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                  padj = replace(padj, padj == 0, .Machine$double.xmin)) %>%
    dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"),
                  regulation = factor(regulation, levels = c("no", "up", "down"))) %>%
    dplyr::arrange(regulation)
  
  # labels data
  plot_df_labels <-
    plot_df %>%
    dplyr::filter(regulation != "no") %>%
    dplyr::mutate(gene_name = ifelse(is.na(gene_name), ID, gene_name))
  
  # annotation table
  annotations <- tibble(xpos = Inf,
                        ypos = -Inf,
                        annotateText = str_c("label cutoff: ",
                                             "p-adjusted <= ", 0.05))
  
  # plot
  ma_plot <-
    ggplot() +
    geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), size = 5, shape = 20) +
    # scale_x_log10(limits = c(0.01, results_limits$x_limit),
    #               breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
    #                    breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
    scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                        values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
    scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  # add labels
  ma_plot_labeled <-
    ma_plot +
    # geom_text(data = plot_df_labels,
    #           aes(x = mean, y = lfc, label = gene_name),
    #           check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
    #           colour = "black", fontface = "plain") +
    geom_text_repel(data = plot_df_labels, 
                    aes(x = mean, y = lfc, label = gene_name),
                    color = "black", fontface = "plain", size = 3,
                    box.padding = 0, point.padding = 0.5, segment.color = "grey50", 
                    min.segment.length = 10) +
    geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText),
              colour = "black", fontface = "italic", size = 2.5,
              hjust = 1.03, vjust = -0.5) +
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = c("gray50", "red2", "#1a75ff"))),
           alpha = F) +
    xlab("log2 mean expression") +
    ylab(str_c("log2 fold change: ", "Cre pos.", " / ", "Cre neg.", " ", result, "\n") %>% str_replace_all(., "_", " ")) +
    theme(axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13, angle = 90)) +
    theme(legend.position = "bottom")
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("Freimer_2017", "MA_plot", "CrePos_vs_CreNeg", result, "png", sep = ".")),
         plot = ma_plot_labeled, width = 10, height = 10)
  
}))


### get 3' UTRs of significant genes
## load mart from Ensembl Biomart
# set animal
ensembl_name <- "mmusculus"

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
                         "http://sep2019.archive.ensembl.org",
                         "http://apr2019.archive.ensembl.org",
                         "http://jan2019.archive.ensembl.org",
                         "http://oct2018.archive.ensembl.org",
                         "http://jul2018.archive.ensembl.org",
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_version) %$%
  URL_archive

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)


### get sequences
# get table of transcripts
utrs_tb <- 
  results_significant_list[["miR-15a"]] %>% 
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::filter(!is.na(gene_id)) %>%
  left_join(., transcripts_info, by = "gene_id") %>% 
  dplyr::select(gene_id, transcript_id, gene_name) %>% 
  dplyr::distinct(.)

# get 3'UTR sequences
UTR_seq <- 
  getSequence(id = utrs_tb$transcript_id, 
              type = "ensembl_transcript_id", 
              seqType = "3utr", 
              mart = mart) %>% 
  as_tibble(.) %>%
  dplyr::left_join(., utrs_tb %>% dplyr::select(transcript_id, gene_name), by = c("ensembl_transcript_id" = "transcript_id")) %>% 
  dplyr::select(gene_name, transcript_id = ensembl_transcript_id, seq_3UTR = `3utr`) %>%
  dplyr::filter(seq_3UTR != "Sequence unavailable") %>%
  dplyr::arrange(gene_name) %>% 
  dplyr::distinct(gene_name, seq_3UTR, .keep_all = T) %>% 
  dplyr::mutate(seq_length = nchar(seq_3UTR)) %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::top_n(., 1, seq_length) %>% 
  dplyr::ungroup(.)

# DNAStringSet
UTR_biostrings <- 
  UTR_seq$seq_3UTR %>% 
  Biostrings::DNAStringSet(.)
names(UTR_biostrings) <- str_c(UTR_seq$gene_name, UTR_seq$transcript_id, sep = ".")

# save fasta
Biostrings::writeXStringSet(x = UTR_biostrings, 
                            filepath = file.path(outpath, 
                                                 str_c("Freimer_2017", "diff_exp", "downregulated", 
                                                       "CrePos_vs_CreNeg", "miR-15a", "3pUTR", "fasta", sep = ".")))

### find miR-15a seed
# seed = 5'-GCTGCT or 5'-TGCTGCT
# get seed
hex_seed <- "GCTGCT"
hepta_seed <- "TGCTGCT"

# find hexamer seed positions
hex_tb <- 
  Biostrings::vmatchPattern(pattern = hex_seed,
                            subject = UTR_biostrings) %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::mutate(hex_position = str_c(start, end, sep = "-")) %>% 
  dplyr::select(gene_name = names, hex_position) %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::summarise(hex_position = str_c(hex_position, collapse = ", "), 
                   hex_count = n())

# find heptamer seed positions
hepta_tb <- 
  Biostrings::vmatchPattern(pattern = hepta_seed,
                            subject = UTR_biostrings) %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::mutate(hepta_position = str_c(start, end, sep = "-")) %>% 
  dplyr::select(gene_name = names, hepta_position) %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::summarise(hepta_position = str_c(hepta_position, collapse = ", "), 
                   hepta_count = n())

# join hexamers and heptamers seeds
seed_tb <- 
  hex_tb %>% 
  dplyr::left_join(., hepta_tb, by = "gene_name") %>% 
  dplyr::select(gene_name, hex_count, hepta_count, hex_position, hepta_position) %>% 
  dplyr::mutate(gene_name = str_remove(gene_name, "\\..*"))

# join with list of differentially expressed genes
results_seed_tb <- 
  results_significant_list[["miR-15a"]] %>% 
  dplyr::filter(log2FoldChange < 0) %>% 
  dplyr::left_join(., seed_tb, by = "gene_name") %>% 
  dplyr::mutate(hex_count = replace(hex_count, is.na(hex_count), 0), 
                hepta_count = replace(hepta_count, is.na(hepta_count), 0)) %>% 
  dplyr::select(ID, gene_name, 
                log2FoldChange, log2_mean_exp, padj, 
                hex_count, hepta_count, hex_position, hepta_position,
                gene_description, coordinates) %>% 
  dplyr::arrange(log2FoldChange) %T>%
  readr::write_csv(., file.path(outpath, 
                                str_c("Freimer_2017", "diff_exp", "significant", "CrePos_vs_CreNeg", "miR_15a.seeds", "csv", sep = ".")))


# ### get all mouse 3' UTRs for Sylamer analysis
# # get table of transcripts
# utrs_tb_all <- 
#   results_list[["miR-15a"]] %>% 
#   dplyr::filter(!is.na(gene_id)) %>% 
#   left_join(., transcripts_info, by = "gene_id") 
# 
# # get 3'UTR sequences
# utrs_seq_all <- getSequence(id = unique(utrs_tb_all$transcript_id), 
#                             type = "ensembl_transcript_id", 
#                             seqType = "3utr", 
#                             mart = mart) 
# 
# # tidy 3'UTR sequences, get the longest 3'UTR per gene
# utrs_tb_tidy <- 
#   utrs_tb_all %>% 
#   left_join(., utrs_seq_all, by = c("transcript_id" = "ensembl_transcript_id")) %>% 
#   dplyr::rename(seq_3UTR = `3utr`) %>% 
#   dplyr::filter(seq_3UTR != "Sequence unavailable") %>%
#   dplyr::mutate(seq_length = nchar(seq_3UTR)) %>% 
#   dplyr::group_by(gene_name) %>% 
#   dplyr::mutate(seq_length_rank = rank(-seq_length, ties.method = "random")) %>%
#   dplyr::filter(seq_length_rank == 1) %>% 
#   dplyr::ungroup(.) %>% 
#   dplyr::select(-seq_length_rank) %>% 
#   dplyr::arrange(-log2FoldChange)
# 
# # DNAStringSet
# UTR_seq_biostrings <- 
#   utrs_tb_tidy$seq_3UTR %>% 
#   Biostrings::DNAStringSet(.)
# names(UTR_seq_biostrings) <- utrs_tb_tidy$gene_name
# 
# # save fasta
# Biostrings::writeXStringSet(x = UTR_seq_biostrings, 
#                             filepath = file.path(outpath,
#                                                  "Sylamer_analysis",
#                                                  str_c("Freimer_2017", "all_probes", "ensembl_93", "3pUTR", "fasta", sep = ".")))
# 
# # save list of gene names ordered by logFC from upregulated to downregulated
# utr_list <- 
#   utrs_tb_tidy %$%
#   gene_name %T>%
#   readr::write_lines(., file.path(outpath, 
#                                   "Sylamer_analysis", 
#                                   "Freimer_2017.all_probes.ensembl_93.3pUTR.logFC_ordered.txt"))
# 
