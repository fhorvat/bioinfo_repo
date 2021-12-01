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

library(limma)
library(Biobase)
library(BiocGenerics)
library(parallel)
library(affyPLM)
library(GEOquery)

library(ggrepel)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### experiment
# set experiment
experiment <- "Joshi_2007_BMCDevBiol_GSE5558"

# base experiment path
base_experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays"

# experiment path 
experiment_path <- file.path(base_experiment_path, experiment)

# set inpath 
inpath <- file.path(experiment_path)

# expressionSet object mad from series matrix from GEO path
exp_set_path <- list.files(inpath, "_series_matrix\\.txt\\.gz", full.names = T)

# chip annotation file path
chip_annotation_path <- file.path(inpath, "GPL3771.annot")

# sample table path
sample_table_path <- list.files(inpath, "\\.sampleTable\\.csv", full.names = T)


### analysis
# base analysis path
base_analysis_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/diffExp"

# analysis path
analysis_path <- file.path(base_analysis_path, experiment)

# set working dir
setwd(analysis_path)

# set outpath
outpath <- getwd()


######################################################## READ DATA
# read expressionSet object mad from series matrix from GEO path
eset <- getGEO(filename = exp_set_path)

# read chip annotation
chip_annotation <- readr::read_delim(chip_annotation_path, delim = "\t", skip = 27)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
### prepare data
# set KO gene, chip and normalization
gene_KO <- "Figla_Null"
chip <- "VMSRMus20K"
normalization <- "quantiles"

# filter table to include only newborn
sample_table %<>%
  dplyr::filter(stage == "newborn") 

# get connection between probe ID's and gene symbols
chip_tidy <- 
  chip_annotation %>% 
  dplyr::select(ID, gene_symbol = `Gene symbol`, gene_name = `Gene title`, coordinates = `Chromosome annotation`) %>% 
  tidyr::separate(coordinates, into = c("coordinates", "coordinates2"), sep = "///") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("coordinates")), 
                   .funs = list(~str_replace(., "Chromosome ", "chr") %>% 
                                  str_replace(",.+ \\(", " ") %>% 
                                  str_replace("\\.\\.", " ") %>% 
                                  str_remove_all("\\)|, complement")))


### normalize expression values
# normalize
eset_norm <-
  eset %>%
  affyPLM::normalize.ExpressionSet.quantiles(.)

# # check normalization
# boxplot(exprs(eset_norm))


#### fit model and get results
# subset expression set
eset_norm <- eset_norm[, sampleNames(eset_norm) %in% sample_table$geo_accession]

# create design matrix
f <- factor(sample_table$genotype, levels = c(gene_KO, "WT"))
design <- model.matrix(~0 + f)
colnames(design) <- c(gene_KO, "WT")
contrast.matrix <- makeContrasts(str_c(gene_KO, " - WT"), levels = design)

# fit linear model
fit <- 
  lmFit(eset_norm, design) %>% 
  contrasts.fit(., contrast.matrix) %>% 
  eBayes(.)


### get and save results
# get the results, tidy
results_tb <- 
  topTable(fit, n = nrow(exprs(eset_norm))) %>% 
  as_tibble(.) %>% 
  dplyr::select(ID, log2FC = logFC, log2_mean_exp = AveExpr, padj = adj.P.Val) %>% 
  dplyr::left_join(., chip_tidy, by = "ID") %>% 
  dplyr::arrange(padj)

# write results to xlsx
wb_all <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb_all, sheetName = str_c(gene_KO, "_vs_WT"))
openxlsx::writeData(wb = wb_all, sheet = str_c(gene_KO, "_vs_WT"), x = results_tb)
openxlsx::saveWorkbook(wb = wb_all, 
                       file = file.path(outpath, str_c("diffExp", str_c(gene_KO, "_vs_WT"), "newborn_ovaries", 
                                                       str_c(normalization, "_normalized"), "limma", chip,
                                                       "all_results.xlsx", sep = ".")), 
                       overwrite = TRUE)

# statistically significant probes
results_tb_sign <-
  results_tb %>%
  dplyr::filter(padj < 0.05)

# check and write
if(nrow(results_tb_sign) > 0){
  
  # write data
  wb_significant <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb_significant, sheetName = str_c(gene_KO, "_vs_WT"))
  openxlsx::writeData(wb = wb_significant, sheet = str_c(gene_KO, "_vs_WT"), x = results_tb_sign)
  openxlsx::saveWorkbook(wb = wb_significant, 
                         file = file.path(outpath, str_c("diffExp", str_c(gene_KO, "_vs_WT"), "newborn_ovaries", 
                                                         str_c(normalization, "_normalized"), "limma", chip,
                                                         "significant_results.xlsx", sep = ".")), 
                         overwrite = TRUE)
  
}


### plot MA plot
# data for plot
plot_df <- 
  results_tb %>% 
  dplyr::select(mean = log2_mean_exp, lfc = log2FC, padj, gene_symbol, ID, gene_name, coordinates) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, padj > 0.05, "no")) %>%
  dplyr::mutate(regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
  dplyr::arrange(regulation)

# get annotation data
annotation_df <- 
  plot_df %>% 
  dplyr::filter(gene_symbol %in% c("C86187"))

# plot
ma_plot <- 
  ggplot() + 
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), shape = 20, size = 5) +
  geom_point(data = annotation_df, aes(x = mean, y = lfc), color = "black", alpha = 1, size = 5, shape = 20) +
  scale_x_continuous(limits = c(-4, 4), breaks = c(-4:4)) +
  scale_y_continuous(limits = c(-12.5, 12.5), breaks = c(-12:12)) +
  scale_colour_manual(values = c(no = "gray50", down = "#1a75ff", up = "red3")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none")

# save plot
ggsave(filename = file.path(outpath, str_c("diffExp", str_c(gene_KO, "_vs_WT"), "newborn_ovaries", 
                                           str_c(normalization, "_normalized"), "limma", chip,
                                           "MA_plot.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)

# add lables to 3 genes
ma_plot <- 
  ma_plot + 
  geom_label_repel(data = annotation_df, aes(x = mean, y = lfc, label = gene_symbol), 
                   fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")

# save plot
ggsave(filename = file.path(outpath, str_c("diffExp", str_c(gene_KO, "_vs_WT"), "newborn_ovaries", 
                                           str_c(normalization, "_normalized"), "limma", chip,
                                           "MA_plot.labels.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)

