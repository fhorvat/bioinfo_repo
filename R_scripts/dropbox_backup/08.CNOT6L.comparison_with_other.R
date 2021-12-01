### INFO: produce scatterplots of ratio of expression in Fugaku and CNOT6L data
### DATE: Sun Mar 11 00:36:34 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/target_comparison")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(VennDiagram)
library(geneplotter)
library(GeneOverlap)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main outpath
outpath <- getwd()

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv"

# FPKM and results path
results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis"
  
# differential expression results
diffExp_results <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/diffExp"

# results Morgan path
Morgan_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/target_comparison"

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

# CNOT6L significantly diff. exp. genes
results_CNOT6L <- 
  list.files(path = diffExp_results, pattern = "diffExp.CNOT6L.*.signif.csv", full.names = T) %>% 
  lapply(., readr::read_csv) %>% 
  dplyr::bind_rows(.)

# CNOT6L FPKM
fpkm_CNOT6L <- readr::read_csv(file = file.path(results_path, "ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"))

# Fugaku FPKM
fpkm_Fugaku <- readr::read_csv(file = file.path(results_path, "ensembl.GRCm38.89.Fugaku.FPKM.csv"))

# Ma 2013 significantly diff. exp. genes
results_Ma <- readr::read_csv(file.path(results_path, "affy.mouse4302.Ma2013.controlVsDcp1a_Dcp2.upregulated.csv"))

# Ivanova 2017 significantly upregulated genes
results_Ivanova <- readr::read_csv(file.path(results_path, "affy.MoGene2.0st.Ivanova2017.KOvsWT.upregulated.csv"))

# Yu 2016 significantly upregulated genes
results_Yu <- 
  list.files(file.path(results_path, "results"), "diffExp.Yu2016.*.signif.csv", full.names = T) %>% 
  lapply(., readr::read_csv) %>% 
  dplyr::bind_rows(.)

# Morgan 2017 significantly upregulated genes
results_Morgan <- readr::read_csv(file.path(Morgan_path, "Morgan2017.affy.MoGene2.0st.KOvsWT.upregulated.csv"))

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_coding <- 
  ensembl_genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# filter results table
results <- 
  results_CNOT6L %>% 
  dplyr::mutate(stage = str_extract(comparison, "1C|MII|GV")) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::select(gene_id, stage) %>% 
  dplyr::mutate(signif = "yes")

# join tables
fpkm_df <- 
  dplyr::left_join(fpkm_CNOT6L, fpkm_Fugaku, by = "gene_id") %>% 
  dplyr::filter(gene_id %in% protein_coding) %>% 
  data.table::setnames(., 2:ncol(.), str_c("s.", colnames(.)[2:ncol(.)]))

# filter Ma
results_Ma %<>%
  dplyr::filter(gene_id %in% protein_coding) %>% 
  dplyr::rename(logFC_Ma = logFC)

# filter Ivanova
results_Ivanova %<>% 
  dplyr::filter(gene_id %in% protein_coding) %>% 
  dplyr::rename(logFC_Ivanova = logFC)

# filter Yu
results_Yu %<>%
  dplyr::filter(gene_id %in% protein_coding, 
                log2FoldChange > 0) %>% 
  dplyr::rename(logFC_Yu = log2FoldChange)

# filter Morgan results
results_Morgan %<>%
  dplyr::filter(gene_id %in% protein_coding) %>% 
  dplyr::rename(logFC_Morgan = logFC)


#### comparison with other KOs (Ma Dcpp, Ivanova YTHDF2, Yu BTG4, Morgan Tut4/7)
### MII / GV with Ma Dcpp KO 
# set significant stage
signif_stage <- "MII"

## MA plot
ma_plot <- 
  fpkm_df %>% 
  dplyr::select(gene_id, s.GV_WT, s.MII_WT) %>% 
  dplyr::mutate(mean = ((s.GV_WT + s.MII_WT) / 2),
                lfc = log2(s.MII_WT + 1) - log2(s.GV_WT + 1)) %>% 
  dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>% 
  dplyr::mutate(signif = replace(signif, is.na(signif), "no"), 
                Ma_upregulated = (gene_id %in% results_Ma$gene_id), 
                CNOT6L_upregulated = (signif == "yes")) %>% 
  mutate_cond(Ma_upregulated, upregulation = "Ma") %>% 
  mutate_cond(CNOT6L_upregulated, upregulation = "CNOT6L") %>% 
  mutate_cond(Ma_upregulated & CNOT6L_upregulated, upregulation = "both") %>% 
  mutate_cond(!(Ma_upregulated | CNOT6L_upregulated), upregulation = "neither") %>% 
  dplyr::mutate(upregulation = factor(upregulation, levels = c("neither", "CNOT6L", "Ma", "both"))) %>%
  dplyr::arrange(upregulation) %>% 
  dplyr::select(-c(s.MII_WT, s.GV_WT, signif, Ma_upregulated, CNOT6L_upregulated)) %>% 
  ggplot(data = ., aes(x = mean, y = lfc, color = upregulation)) + 
  geom_point(size = 5, shape = 20) +
  scale_x_log10(limits = c(0.001, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(breaks = c(-2, 2)) +
  scale_color_manual(values = c(neither = "gray60", CNOT6L = "red3", Ma = "black", both = "#1a75ff")) +
  guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2(MII WT / GV WT)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
        axis.title.y = element_text(size = 15, vjust = 0.3), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot.FPKM.CNOT6L.MII_vs_GV.WT.Ma.", signif_stage, ".KO.sign.GRCm38.89.png")),
       plot = ma_plot, width = 15, height = 10)

## Venn diagram
png(file = file.path(outpath, "Venn.Ma_vs_CNOT6L.significant.upregulated.png"),
    width = 1000, height = 1000, units = "px", type = "cairo")
up_plot <- venn.diagram(x = list(Ma = na.omit(results_Ma$gene_id),
                                 CNOT6L = results %>% dplyr::filter(stage == signif_stage, signif == "yes") %$% gene_id),
                        filename = NULL,
                        fill = c("blue", "red"),
                        alpha = c(0.5, 0.5),
                        cex = 4,
                        cat.cex = 2,
                        main = "MII KO vs. WT upregulated genes",
                        main.cex = 2)
grid.draw(up_plot)
dev.off()


### MII / GV with Ivanova YTHDF2 KO
# set significant stage
signif_stage <- "MII"

## MA plot
ma_plot <- 
  fpkm_df %>% 
  dplyr::select(gene_id, s.GV_WT, s.MII_WT) %>% 
  dplyr::mutate(mean = ((s.GV_WT + s.MII_WT) / 2),
                lfc = log2(s.MII_WT + 1) - log2(s.GV_WT + 1)) %>% 
  dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>% 
  dplyr::mutate(signif = replace(signif, is.na(signif), "no"), 
                Ivanova_upregulated = (gene_id %in% results_Ivanova$gene_id), 
                CNOT6L_upregulated = (signif == "yes")) %>% 
  mutate_cond(Ivanova_upregulated, upregulation = "Ivanova") %>% 
  mutate_cond(CNOT6L_upregulated, upregulation = "CNOT6L") %>% 
  mutate_cond(Ivanova_upregulated & CNOT6L_upregulated, upregulation = "both") %>% 
  mutate_cond(!(Ivanova_upregulated | CNOT6L_upregulated), upregulation = "neither") %>% 
  dplyr::mutate(upregulation = factor(upregulation, levels = c("neither", "CNOT6L", "Ivanova", "both"))) %>%
  dplyr::arrange(upregulation) %>% 
  dplyr::select(-c(s.MII_WT, s.GV_WT, signif, Ivanova_upregulated, CNOT6L_upregulated)) %>% 
  ggplot(data = ., aes(x = mean, y = lfc, color = upregulation)) + 
  geom_point(size = 5, shape = 20) +
  scale_x_log10(limits = c(0.001, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(breaks = c(-2, 2)) +
  scale_color_manual(values = c(neither = "gray60", CNOT6L = "red3", Ivanova = "black", both = "#1a75ff")) +
  guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2(MII WT / GV WT)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
        axis.title.y = element_text(size = 15, vjust = 0.3), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot.FPKM.CNOT6L.MII_vs_GV.WT.Ivanova.", signif_stage, ".KO.sign.GRCm38.89.png")),
       plot = ma_plot, width = 15, height = 10)

## Venn diagram
png(file = file.path(outpath, "Venn.Ivanova_vs_CNOT6L.significant.upregulated.png"),
    width = 1000, height = 1000, units = "px", type = "cairo")
up_plot <- venn.diagram(x = list(YTHDF2 = na.omit(results_Ivanova$gene_id),
                                 CNOT6L = results %>% dplyr::filter(stage == signif_stage, signif == "yes") %$% gene_id),
                        filename = NULL,
                        fill = c("blue", "red"),
                        alpha = c(0.5, 0.5),
                        cex = 4,
                        cat.cex = 2,
                        main = "MII KO vs. WT upregulated genes",
                        main.cex = 2)
grid.draw(up_plot)
dev.off()


### MII / GV with Yu BTG4 KO
# set significant stage
signif_stage <- "MII"

## MA plot
ma_plot <- 
  fpkm_df %>% 
  dplyr::select(gene_id, s.GV_WT, s.MII_WT) %>% 
  dplyr::mutate(mean = ((s.GV_WT + s.MII_WT) / 2),
                lfc = log2(s.MII_WT + 1) - log2(s.GV_WT + 1)) %>% 
  dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>% 
  dplyr::mutate(signif = replace(signif, is.na(signif), "no"), 
                Yu_upregulated = (gene_id %in% (results_Yu %>% dplyr::filter(str_detect(.$comparison, signif_stage)) %$% gene_id)), 
                CNOT6L_upregulated = (signif == "yes")) %>% 
  mutate_cond(Yu_upregulated, upregulation = "Yu") %>% 
  mutate_cond(CNOT6L_upregulated, upregulation = "CNOT6L") %>% 
  mutate_cond(Yu_upregulated & CNOT6L_upregulated, upregulation = "both") %>% 
  mutate_cond(!(Yu_upregulated | CNOT6L_upregulated), upregulation = "neither") %>% 
  dplyr::mutate(upregulation = factor(upregulation, levels = c("neither", "Yu", "CNOT6L", "both"))) %>%
  dplyr::arrange(upregulation) %>% 
  dplyr::select(-c(s.MII_WT, s.GV_WT, signif, Yu_upregulated, CNOT6L_upregulated)) %>% 
  ggplot(data = ., aes(x = mean, y = lfc, color = upregulation)) + 
  geom_point(size = 5, shape = 20) +
  scale_x_log10(limits = c(0.001, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(breaks = c(-2, 2)) +
  scale_color_manual(values = c(neither = "gray60", CNOT6L = "red3", Yu = "black", both = "#1a75ff")) +
  guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2(MII WT / GV WT)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
        axis.title.y = element_text(size = 15, vjust = 0.3), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot.FPKM.CNOT6L.MII_vs_GV.WT.Yu.", signif_stage, ".KO.sign.GRCm38.89.png")),
       plot = ma_plot, width = 15, height = 10)


## Venn diagram
png(file = file.path(outpath, str_c("Venn.Yu_vs_CNOT6L.", signif_stage, ".significant.upregulated.png")),
    width = 1000, height = 1000, units = "px", type = "cairo")
up_plot <- venn.diagram(x = list(BTG4 = results_Yu %>% dplyr::filter(str_detect(comparison, signif_stage)) %$% gene_id,
                                 CNOT6L = results %>% dplyr::filter(stage == signif_stage) %$% gene_id),
                        filename = NULL,
                        fill = c("blue", "red"),
                        alpha = c(0.5, 0.5),
                        cex = 4,
                        cat.cex = 2,
                        main = str_c(signif_stage, " BTG4/CNOT6L KO vs. WT upregulated genes"),
                        main.cex = 2)
grid.draw(up_plot)
dev.off()


### MII / GV with Morgan Tut4/7 KO
# set significant stage
signif_stage <- "MII"

## MA plot
ma_plot <- 
  fpkm_df %>% 
  dplyr::select(gene_id, s.GV_WT, s.MII_WT) %>% 
  dplyr::mutate(mean = ((s.GV_WT + s.MII_WT) / 2),
                lfc = log2(s.MII_WT + 1) - log2(s.GV_WT + 1)) %>% 
  dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>% 
  dplyr::mutate(signif = replace(signif, is.na(signif), "no"), 
                Morgan_upregulated = (gene_id %in% results_Morgan$gene_id), 
                CNOT6L_upregulated = (signif == "yes")) %>% 
  mutate_cond(Morgan_upregulated, upregulation = "Morgan") %>% 
  mutate_cond(CNOT6L_upregulated, upregulation = "CNOT6L") %>% 
  mutate_cond(Morgan_upregulated & CNOT6L_upregulated, upregulation = "both") %>% 
  mutate_cond(!(Morgan_upregulated | CNOT6L_upregulated), upregulation = "neither") %>% 
  dplyr::mutate(upregulation = factor(upregulation, levels = c("neither", "CNOT6L", "Morgan", "both"))) %>%
  dplyr::arrange(upregulation) %>% 
  dplyr::select(-c(s.MII_WT, s.GV_WT, signif, Morgan_upregulated, CNOT6L_upregulated)) %>% 
  ggplot(data = ., aes(x = mean, y = lfc, color = upregulation)) + 
  geom_point(size = 5, shape = 20) +
  scale_x_log10(limits = c(0.001, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(breaks = c(-2, 2)) +
  scale_color_manual(values = c(neither = "gray60", CNOT6L = "red3", Morgan = "black", both = "#1a75ff")) +
  guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2(MII WT / GV WT)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
        axis.title.y = element_text(size = 15, vjust = 0.3), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot.FPKM.CNOT6L.MII_vs_GV.WT.MorganGV.", signif_stage, ".KO.sign.GRCm38.89.png")),
       plot = ma_plot, width = 15, height = 10)

## Venn diagram
png(file = file.path(outpath, str_c("Venn.Morgan_vs_CNOT6L.", signif_stage, ".significant.upregulated.png")),
    width = 1000, height = 1000, units = "px", type = "cairo")
up_plot <- venn.diagram(x = list(Morgan = na.omit(results_Morgan$gene_id),
                                 CNOT6L = results %>% dplyr::filter(stage == signif_stage, signif == "yes") %$% gene_id),
                        filename = NULL,
                        fill = c("blue", "red"),
                        alpha = c(0.5, 0.5),
                        cex = 4,
                        cat.cex = 2,
                        main = "CNOT6L MII vs. Morgan GV - KO vs. WT upregulated genes",
                        main.cex = 2)
grid.draw(up_plot)
dev.off()



### statistical test
# test: Fisher exact test
# null hypothesis: there is no association between X and Y (CNOT6L and other KOs are not related)
# If the null hypothesis were true (if CNOT6L and other KOs are not related) how likely is it that 
# overlap between upregulated genes between CNOT6L KO in MII and other KOs (DCPP, YTHDF2, BTG4, Tut4/7) is 
# as large as observed or larger?

## GeneOverlap
newGOM(gsetA = list(CNOT6L = results %>% dplyr::filter(stage == "MII", signif == "yes") %$% gene_id), 
       gsetB = list(DCPP = results_Ma$gene_id, 
                    YTHDF2 = results_Ivanova$gene_id, 
                    BTG4 = results_Yu %>% dplyr::filter(str_detect(comparison, "MII")) %$% gene_id),
       genome.size = length(protein_coding)) %>% 
  getMatrix(., name = "pval") %>% 
  as.data.frame(.) %>% 
  set_colnames(., "p.value") %>% 
  tibble::rownames_to_column(., var = "KO") %>% 
  dplyr::mutate(is.significant = ifelse(p.value <= 0.05, "yes", "no"), 
                KO = stringr::str_remove(KO, ".CNOT6L"))

## base
# set data
CNOT6L_genes <- results %>% dplyr::filter(stage == "MII", signif == "yes") %$% gene_id
DCPP_genes <- results_Ma$gene_id
YTHDF2_genes <- results_Ivanova$gene_id
BTG4_genes <- results_Yu %>% dplyr::filter(str_detect(comparison, "MII")) %$% gene_id
TUT_genes <- results_Morgan$gene_id
all_genes <- protein_coding

# set variables
in_A <- CNOT6L_genes
in_B <- TUT_genes

## create contingency table
# in none: n - union(A, B)
# only in A: setdiff(A, B)
# only in B: setdiff(B, A)
# in both: intersect(A, B)
# matrix(in none, only in A, only in B, in both), nrow = 2)
cont_table <- matrix(c(length(all_genes) - length(union(in_A, in_B)), 
                       length(setdiff(in_A, in_B)), 
                       length(setdiff(in_B, in_A)), 
                       length(intersect(in_A, in_B))), 
                     nrow = 2)

# test
fisher.test(cont_table, alternative = "l")
