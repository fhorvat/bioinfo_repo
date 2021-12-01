### INFO: 
### DATE: Fri Nov 29 12:17:12 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Wu_2019_unpub_GSE133748/Analysis/ERCC_spike_expression")

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

library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# fpkm path
fpkm_path <- file.path(inpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.ERCC.FPKM.csv")

# ERCC info path
ercc_path <- "/common/DB/genome_reference/spike_in/ERCC92/ERCC92.info.txt"

######################################################## READ DATA
# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

# read ERCC info
ercc_tb <- readr::read_delim(ercc_path, delim = "\t")

######################################################## MAIN CODE
# clean ERCC info and calculate expected number of molecules
ercc_molecules <- 
  ercc_tb %>% 
  set_colnames(., c("num", "id", "subgroup", "conc_mix1", "conc_mix2", "expected_fc", "log2_mix1_mix2")) %>% # concentraction is in attomoles/ul
  dplyr::select(id, conc_mix1) %>% 
  dplyr::mutate(ercc_molecules_n = 
                  (conc_mix1 / 100000) * # mix was diluted 10^5 fold
                  1 *                    # 1 µl of the diluted ERCC mix was added to each sample
                  1/10^18 *              # Number of attomoles in a mole
                  6.02214179e23)         # Number of molecules in a mole

# join with ERCC FPKM
ercc_fpkm <- 
  ercc_molecules %>% 
  dplyr::select(-conc_mix1) %>%
  dplyr::left_join(., fpkm_tb %>% dplyr::select_at(vars(gene_id, starts_with("s_"))), by = c("id" = "gene_id")) %>% 
  tidyr::pivot_longer(-c(id, ercc_molecules_n),
                      names_to = "sample_id", 
                      values_to = "fpkm") %>% 
  dplyr::filter(ercc_molecules_n >= 1, 
                fpkm > 0) %>% 
  dplyr::select(gene_name = id, sample_id, fpkm, molecules_n = ercc_molecules_n)

### linear model
# build linear regression model - line has to go through origin so the zero term is added to the formula
linearMod <- lm(molecules_n ~ 0 + fpkm, data = ercc_fpkm)  # build linear regression model on full data
summary(linearMod)

# calculate predictions intervals
pred.int <- predict(linearMod, interval = "prediction")

# add to the table
ercc_fpkm <- 
  cbind(ercc_fpkm, pred.int) %>% 
  as_tibble(.)

# get annotation table
annotations <- tibble(xpos = Inf,
                      ypos = -Inf,
                      annotateText = str_c("R^2 adjusted ", round(summary(linearMod)$adj.r.squared, 4)))

# plot regression line + confidence intervals
ln_plot <- 
  ggplot() +
  geom_point(data = ercc_fpkm, aes(x = fpkm, y = molecules_n)) +
  stat_smooth(data = ercc_fpkm, aes(x = fpkm, y = molecules_n), method = lm) +
  geom_line(data = ercc_fpkm, aes(x = fpkm, y = lwr), color = "red", linetype = "dashed")+
  geom_line(data = ercc_fpkm, aes(x = fpkm, y = upr), color = "red", linetype = "dashed") +
  geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText),
            colour = "black", fontface = "italic", size = 2.5,
            hjust = 1.03, vjust = -0.5) +
  coord_cartesian(xlim = c(0, max(ercc_fpkm$fpkm)),
                  ylim = c(0, max(ercc_fpkm$molecules_n))) +
  xlab("FPKM") +
  ylab("ERCC no. molecules") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  ggsave(filename = file.path(outpath, "ERCC.n_molecules.fpkm.linear_regression.png"),
         height = 10, width = 12)


### get estimated number of molecules for all genes
# tidy FPKM table
count_tb <- 
  fpkm_tb %>% 
  dplyr::select(gene_id, starts_with("s_")) %>% 
  tidyr::pivot_longer(-c(gene_id),
                      names_to = "sample_id", 
                      values_to = "fpkm") %>% 
  dplyr::mutate(molecules_n = predict(linearMod, newdata = .), 
                molecules_n = round(molecules_n, 0)) %>%
  tidyr::pivot_wider(id_cols = gene_id, names_from = sample_id, values_from = molecules_n) %T>%
  readr::write_csv(., file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.ERCC.Wu_2019.n_molecules.csv"))


### get FPKM of chosen genes, calculate number of molecules, join with ERCC, plot together 
# set chosen gene names
chosen_genes <- 
  c("Drosha" , "Dgcr8" , "Dicer1" , "Tarbp2" , "Ago1" , "Ago2" , "Ago3" , "Ago4" , 
    "Tnrc6a" , "Tnrc6b" , "Tnrc6c" , "Ddx6" , "Nr6a1" , "Hif1an" , "Hmga2" , "Ccnb1" , 
    "Zp3" , "Rps17" , "Mos" , "Plat" , "Dcp1a" , "Rps17" , "Eif3l" , "Ppil3" , 
    "Hprt" , "Kif2c" , "Elob" , "Elobl" , "Dcp1a" , "C86187") %>% 
  unique(.)

# filter FPKM dataset, calculate number of molecules based on the linear regression model of ERCC spike-in
fpkm_tidy <- 
  fpkm_tb %>% 
  dplyr::filter(gene_name %in% chosen_genes) %>% 
  dplyr::mutate(gene_name = replace(gene_name, gene_name == "C86187", "Sirena1")) %>% 
  dplyr::select(gene_name, starts_with("s_")) %>% 
  tidyr::pivot_longer(-c(gene_name),
                      names_to = "sample_id", 
                      values_to = "fpkm") %>% 
  dplyr::mutate(molecules_n = predict(linearMod, newdata = .))

# get mean FPKM and number of molecules in all samples
fpkm_mean <- 
  fpkm_tidy %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::summarise(mean_fpkm = mean(fpkm) %>% round(., 3), 
                   mean_molecules_n = mean(molecules_n) %>% round(., 3)) %>% 
  dplyr::arrange(match(gene_name, replace(chosen_genes, chosen_genes == "C86187", "Sirena1"))) %T>%
  readr::write_csv(., file.path(outpath, "chosen_genes.Wu_2019_unpub_GSE133748.number_of_molecules_from_ERCC.mean.csv"))


### plot with ERCC
# plot
plot_fpkm <- 
  ggplot() +
  geom_point(data = ercc_fpkm, aes(x = fpkm, y = molecules_n), color = "black", size = 1) +
  geom_line(data = ercc_fpkm, aes(x = fpkm, y = fit), color = "blue") +
  geom_point(data = fpkm_mean, aes(x = mean_fpkm, y = mean_molecules_n), color = "red", size = 3, shape = "triangle") +
  geom_label_repel(data = fpkm_mean, aes(x = mean_fpkm, y = mean_molecules_n, label = gene_name)) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("FPKM (log10)") +
  ylab("ERCC predicted no. molecules (log10)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  ggsave(filename = file.path(outpath, "chosen_genes.Wu_2019_unpub_GSE133748.number_of_molecules_from_ERCC.mean.png"),
         height = 10, width = 12)


