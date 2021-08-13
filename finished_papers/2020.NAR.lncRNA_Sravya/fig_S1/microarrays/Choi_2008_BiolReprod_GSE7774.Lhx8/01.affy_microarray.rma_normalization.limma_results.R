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
library(affy)
library(mouse4302.db)

library(ggrepel)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### experiment (chip = Affymetrix Mouse Genome 430 2.0 Array, GPL1261)
# set experiment
experiment <- "Choi_2008_BiolReprod_GSE7774"

# base experiment path
base_experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays"

# experiment path 
experiment_path <- file.path(base_experiment_path, experiment)

# set inpath 
inpath <- file.path(experiment_path)

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
# read sample table
sample_table <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
### prepare data
# set KO gene, chip and normalization
gene_KO <- "LHX8_Null"
chip <- "Mouse430_2"
normalization <- "rma"

# prepare for sample table for AnnotatedDataFrame
sample_table %<>% 
  as.data.frame(.) %>% 
  column_to_rownames(var = "sample_id")

# set sample names
sample_names <- 
  sample_table$geo_accession %>% 
  magrittr::set_names(rownames(sample_table))

# read affymetrix microarray set
eset <- affy::read.affybatch(filenames = sample_table$cel_path,
                             phenoData = AnnotatedDataFrame(sample_table))

# normalize using RMA algorithm
eset_norm <- rma(eset)


#### fit model and get results
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
  as_tibble(., rownames = "ID") %>% 
  dplyr::select(ID, log2FC = logFC, log2_mean_exp = AveExpr, padj = adj.P.Val) %>% 
  dplyr::mutate(gene_id = mapIds(mouse4302.db, keys = ID, column = "ENSEMBL", keytype = "PROBEID", multiVals = "first"), 
                gene_symbol = mapIds(mouse4302.db, keys = ID, column = "SYMBOL", keytype = "PROBEID", multiVals = "first"), 
                gene_name = mapIds(mouse4302.db, keys = ID, column = "GENENAME", keytype = "PROBEID", multiVals = "first")) %>% 
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


### plot MA plot
# data for plot
plot_df <- 
  results_tb %>% 
  dplyr::select(mean = log2_mean_exp, lfc = log2FC, padj, gene_symbol, ID, gene_name) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, padj > 0.05, "no")) %>%
  dplyr::mutate(regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
  dplyr::arrange(regulation)

# get annotation data
annotation_df <- 
  plot_df %>% 
  dplyr::filter(gene_symbol %in% c("C86187", "Nobox", "Lhx8", "Figla"))

# plot
ma_plot <- 
  ggplot() + 
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), shape = 20, size = 5) +
  geom_point(data = annotation_df, aes(x = mean, y = lfc), color = "black", alpha = 1, size = 5, shape = 20) +
  scale_x_continuous(limits = c(2, 15), breaks = c(2:15)) +
  scale_y_continuous(limits = c(-8, 8), breaks = c(-8:8)) +
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

