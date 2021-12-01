### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### experiment
# set experiment
experiment <- "Freimer_2017_CurrBiol_GSE92658"

# base experiment path
base_experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays"

# experiment path 
experiment_path <- file.path(base_experiment_path, experiment)

# set working dir
setwd(experiment_path)

# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


## chip files
# list normalized .csv files
norm_files_list <- list.files(path = inpath, pattern = "norm\\.csv\\.gz", full.names = T)

# chip annotation file path
chip_annotation_path <- file.path(inpath, "GPL6885.annot.gz")


## gene annoation
# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

######################################################## READ DATA
# read .csv files
norm_tb <- purrr::map(norm_files_list, function(path){
  
  # read file
  tb <- 
    readr::read_csv(path) %>% 
    dplyr::select(-contains("Pval")) %>% 
    dplyr::rename(probe = X1)
  
}) %>% 
  purrr::reduce(., left_join, by = "probe") %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "probe") %>% 
  as.matrix(.)

# read chip annotation
chip_annotation <- readr::read_delim(chip_annotation_path, delim = "\t", skip = 28)

# read genes info
genes_info <- 
  readr::read_csv(genes_info_path) %>% 
  tidyr::unite(., coordinates, seqnames, start, end, sep = " ")

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
                                  str_remove_all("\\)|, complement")))

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
                  replace(., . == "none", "ctrl") %>% 
                  str_replace(., "miR-15a", "miR15")) %>% 
  dplyr::select(sample_id, geo_accession, stage, genotype, treatment)


### get differentially expressed genes
# get normalized expression, change column order to match sample table
exprs_tb <- norm_tb[, match(sample_table$sample_id, colnames(norm_tb))]

# check
if(all(colnames(exprs_tb) != sample_table$sample_id)){
  stop("Order of columns doesn't match order of samples in sample table")
}

# create design matrix
TS <- paste(sample_table$genotype, sample_table$treatment, sep = ".")
TS <- factor(TS, levels = c("Cre_Neg.ctrl", "Cre_Pos.ctrl", "Cre_Neg.miR15", "Cre_Pos.miR15"))
design <- model.matrix(~0 + TS)
colnames(design) <- levels(TS)

# fit the model
fit <- lmFit(exprs_tb, design)
cont.matrix <- makeContrasts(Cre_PosvsCreNeg_in_ctrl = Cre_Pos.ctrl - Cre_Neg.ctrl, 
                             Cre_PosvsCreNeg_in_miR15a = Cre_Neg.miR15 - Cre_Pos.miR15, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

results <- decideTests(fit2)

# get the results, tidy
results_tb <- 
  topTable(fit2, n = nrow(exprs_tb), coef = "Cre_PosvsCreNeg_in_miR15a") %>% 
  as_tibble(., rownames = "ID") %>% 
  dplyr::select(ID, log2FoldChange = logFC, log2_mean_exp = AveExpr, padj = adj.P.Val) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::left_join(., chip_tidy %>% dplyr::select(ID, gene_name), by = "ID") %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_name, gene_description, coordinates), by = "gene_name")
# readr::write_csv(., file.path(outpath, str_c("Freimer_2017", "diff_exp", "CrePos_vs_CreNeg", result, "csv", sep = ".")))

# filter table
results_df_sign <-
  results_tb  %>%
  dplyr::filter(padj <= 0.05)


### plot MA plots
## prepare results
# get results table
results_df <- results_tb

# significant results
results_df_sign <-
  results_df_sign %>%
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down")) %>%
  dplyr::select(ID, regulation)


## MA plot
# set result
result <- "miR15"
  
# data for plot
plot_df <-
  results_df %>%
  dplyr::select(mean = log2_mean_exp, lfc = log2FoldChange, padj, ID, gene_name) %>%
  dplyr::mutate(lfc = -lfc) %>% 
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
  geom_text(data = plot_df_labels,
            aes(x = mean, y = lfc, label = gene_name),
            check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
            colour = "black", fontface = "plain") +
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





