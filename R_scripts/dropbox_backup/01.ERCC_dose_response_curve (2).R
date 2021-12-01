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
# calculate standard error of the mean FPKM value
standard.error <- function(x) {
  sqrt(var(x) / length(x))
}

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

# build linear regression model
linearMod <- lm(molecules_n ~ fpkm, data = ercc_fpkm)  # build linear regression model on full data
summary(linearMod)

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
  dplyr::mutate(molecules_n = round(linearMod$coefficients[1] + (linearMod$coefficients[2] * fpkm)))
dplyr::mutate(molecules_n = (linearMod$coefficients[2] * fpkm))

## get and save wide tables
# FPKM
fpkm_tidy %>%
  pivot_wider(-molecules_n,
              names_from = "sample_id",
              values_from = "fpkm") %T>%
  readr::write_csv(., file.path(outpath, "chosen_genes.Wu_2019_unpub_GSE133748.FPKM.csv"))

# number of molecules
fpkm_tidy %>%
  pivot_wider(-fpkm,
              names_from = "sample_id",
              values_from = "molecules_n") %T>%
  readr::write_csv(., file.path(outpath, "chosen_genes.Wu_2019_unpub_GSE133748.number_of_molecules.from_ERCC.csv"))

# long mean FPKM, number of molecules
fpkm_mean <- 
  fpkm_tidy %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::summarise(mean_fpkm = mean(fpkm), 
                   mean_molecules_n = mean(molecules_n)) %>% 
  dplyr::arrange(match(gene_name, replace(chosen_genes, chosen_genes == "C86187", "Sirena1"))) %T>%
  readr::write_csv(., file.path(outpath, "chosen_genes.Wu_2019_unpub_GSE133748.number_of_molecules.from_ERCC.mean.csv"))


### plot
# calculate mean value of gene FPKM/molecules number
plot_fpkm_tb <- 
  fpkm_mean %>% 
  dplyr::mutate(mean_fpkm = log10(mean_fpkm), 
                mean_molecules_n = log10(mean_molecules_n))

# join with ERCC fpkm, calculate log2 of FPKM and molecule numbers, prepare for plotting
plot_ercc_tb <- 
  ercc_fpkm %>% 
  dplyr::mutate(log2_fpkm = log10(fpkm), 
                log2_molecules_n = log10(molecules_n))

# plot without gene labels
ln_plot <- 
  ggplot() +
  geom_point(data = plot_ercc_tb, 
             mapping = aes(x = log2_molecules_n, y = log2_fpkm), color = "black", size = 3) +
  geom_smooth(data = plot_ercc_tb,
              mapping = aes(x = log2_molecules_n, y = log2_fpkm),
              method = lm) +
  scale_x_continuous(limits = c(0, 6), breaks = 0:6) +
  scale_y_continuous(limits = c(-1.5, 3.5)) +
  xlab("log10(transcript number)") +
  ylab("log10(FPKM)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none")

# save
ggsave(filename = file.path(outpath, "chosen_genes.ERCC.molecules_n.fpkm.log10.plot.png"), plot = ln_plot, height = 10, width = 10)

# add gene labels and plot
ln_plot_labels <- 
  ln_plot +
  geom_label_repel(data = plot_fpkm_tb, 
                   mapping = aes(x = mean_molecules_n, y = mean_fpkm, label = gene_name)) +
  geom_point(data = plot_fpkm_tb, 
             mapping = aes(x = mean_molecules_n, y = mean_fpkm), color = "red3", size = 5)

# save
ggsave(filename = file.path(outpath, "chosen_genes.ERCC.molecules_n.fpkm.log10.plot.labels.png"), plot = ln_plot_labels, height = 10, width = 10)


