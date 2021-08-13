### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/review/germ_cells_count")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# raw counts table path
raw_counts_path <- file.path(inpath, "germ_cell_per_semineferous_tubule.counts.xlsx")

######################################################## READ DATA
# read raw counts data
raw_counts <- openxlsx::read.xlsx(raw_counts_path)

######################################################## MAIN CODE
# tidy table
counts_tb <- 
  raw_counts %>% 
  as_tibble(.) %>% 
  dplyr::rename(genotype = counts) %>% 
  tidyr::pivot_longer(cols = -genotype, names_to = "germcell_present", values_to = "count") %>% 
  dplyr::mutate(germcell_present = replace(germcell_present, germcell_present == "0", "not_present"), 
                germcell_present = replace(germcell_present, germcell_present != "not_present", "present")) %>% 
  dplyr::group_by(genotype, germcell_present) %>% 
  dplyr::summarise(count = sum(count)) %>%
  tidyr::pivot_wider(id_cols = genotype, names_from = germcell_present, values_from = count) %>% 
  as.data.frame(.) %>%
  tibble::column_to_rownames(., var = "genotype")

# chisq test
chisq_result <- chisq.test(x = counts_tb, simulate.p.value = T)

### fit a GLM - Poisson 
# tidy table
counts_tb <- 
  raw_counts %>% 
  as_tibble(.) %>% 
  dplyr::rename(genotype = counts) %>% 
  tidyr::pivot_longer(cols = -genotype, names_to = "germcell_present", values_to = "count") %>% 
  tidyr::uncount(count) %>% 
  dplyr::select(germcell_present, genotype) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "KO")), 
                germcell_present = as.numeric(germcell_present))

# # plot as freqpoly
# plot_bar <- 
#   ggplot() +
#   geom_freqpoly(data = counts_tb, 
#                 mapping = aes(germcell_present, color = genotype), 
#                 binwidth = 1) + 
#   guides(color = guide_legend(reverse = TRUE)) +
#   scale_color_manual(values = c("WT" = "#7f7f7f", "KO" = "#ff0000")) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom")
# 
# # save
# ggsave(plot = plot_bar, 
#        filename = file.path(outpath, "freqpoly.germ_cell_counts.genotype.png"), width = 10, height = 10)

# check mean and variance
mean(counts_tb$germcell_present)
var(counts_tb$germcell_present)

# fit the model
poisson.model <- glm(germcell_present ~ genotype, counts_tb, family = poisson(link = "log"))
summary(poisson.model)

# non-significant for genotype!

# Since Residual deviance is greater than degrees of freedom ver-dispersion exists. 
# This means that the estimates are correct, but the standard errors (standard deviation) 
# are wrong and unaccounted for by the model.
# So, to have a more correct standard error we can use a quasi-poisson model:
poisson.model2 <- glm(germcell_present ~ genotype, counts_tb, family = quasipoisson(link = "log"))
summary(poisson.model2)

### compare the models
library(arm)

# extract coefficients from first model using 'coef()'
coef1 <- coef(poisson.model)

# extract coefficients from second model
coef2 <- coef(poisson.model2)

# extract standard errors from first model using 'se.coef()'
se.coef1 <- se.coef(poisson.model)

# extract standard errors from second model
se.coef2 <- se.coef(poisson.model2)

# use 'cbind()' to combine values into one dataframe
models.both <- cbind(coef1, se.coef1, coef2, se.coef2, exponent = exp(coef1))

# show dataframe
models.both

# In above output, we can see the coefficients are the same, but the standard errors are different.
# Estimate for genotype => value is -0.1205202 and exponent is 0.8864592
# This means that switching genotype from WT to KO results in decrease in number of germ cells in
# tubules 0.8864592 times the intercept. 
# Another way of saying this is if we change genotype from WT to KO, the number of of germ cells in
# tubules will fall by 11.4% (1-0.8864592) assuming all other variables are the same.
 

### Visualizing Findings Using jtools
library(jtools)
library(interactions)

# plot regression coefficients for poisson.model2
png(filename = file.path(outpath, str_c("germ_cell_counts.genotype", "reg_coef", "png", sep = ".")),
    width = 1000, height = 1000)
plot_summs(poisson.model2, scale = TRUE, exp = TRUE)
dev.off()

# plot regression coefficients for poisson.model2 and poisson.model
png(filename = file.path(outpath, str_c("germ_cell_counts.genotype", "reg_coef.2", "png", sep = ".")),
    width = 1000, height = 1000)
plot_summs(poisson.model, poisson.model2, scale = TRUE, exp = TRUE)
dev.off()

# plot interaction between predictor variables (genotype)
png(filename = file.path(outpath, str_c("germ_cell_counts.genotype", "cat_plot", "png", sep = ".")),
    width = 1000, height = 1000)
cat_plot(poisson.model2, pred = genotype, plot.points = TRUE)
dev.off()

