### INFO: 
### DATE: Mon Oct 07 14:20:49 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/rodent_oocytes.small_RNAseq.2021_Sep/Analysis/rat/read_width_classification")

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

library(ggridges)
library(cowplot)
library(scales)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# makes scales nicer (https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot)
fancy_scientific <- function(l) {
  
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  
  # write 0 properly
  l <- gsub("0e\\+00", "0", l)
  
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
  l <- gsub("e\\+", "e", l)
  
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  
  # return this as an expression
  parse(text = l)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/rodent_oocytes.small_RNAseq.2021_Sep"

# mapped path 
mapped_path <- file.path(base_path, "/Data/Mapped/STAR_rn6.deduped")

# list library histogram files
library_hist_path <- list.files(mapped_path, pattern = "library_hist\\..*\\.SE\\.txt", full.names = T, recursive = T)

# list library size files
library_size_path <- list.files(mapped_path, pattern = "library_sizes.txt", full.names = T, recursive = T)

# sample table path
sample_table_path <- file.path(base_path, "/Data/Documentation")
sample_table_path <- list.files(sample_table_path, pattern = ".*\\.sampleTable.csv", full.names = T)

######################################################## READ DATA
# read library histograms 
library_hist_tb <- purrr::map(library_hist_path, function(path){
  
  # read 
  readr::read_delim(file = path, delim = "\t", col_names = c("count", "cigar")) %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% str_remove_all(., "^library_hist\\.|\\.SE\\.perfect\\.txt$")) %>% 
    dplyr::select(sample_id, cigar, count)
  
}) %>% 
  dplyr::bind_rows(.)

# read library size
library_tb <- readr::read_delim(library_size_path, col_names = c("sample_id", "library_size"), delim = "\t")

# read sample table
sample_tb <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
### prepare data
# clean library histogram
library_hist_tidy <- 
  library_hist_tb %>% 
  dplyr::filter(!str_detect(cigar, "I|D")) %>% 
  dplyr::mutate(alignment_length = str_extract(cigar, "[0-9]{2}M") %>% str_remove(., "M"), 
                left_soft = str_extract(cigar, "^[0-9]+S(?=[0-9]{2}M)") %>% str_remove(., "S"), 
                right_soft = str_extract(cigar, "(?<=[0-9]{2}M)[0-9]+S$") %>% str_remove(., "S"), 
                sample_id = str_remove(sample_id, "\\.txt$")) 

# clean library size
library_tb_tidy <- 
  library_tb %>% 
  dplyr::filter(str_detect(sample_id, "\\.19to32nt$")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$")) 


### data for plots
# get read length distribution
read_length_tb <- 
  library_hist_tidy %>% 
  dplyr::select(sample_id, alignment_length, count) %>% 
  dplyr::mutate(alignment_length = as.numeric(alignment_length)) %>% 
  dplyr::group_by(sample_id, alignment_length) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_tb_tidy, by = "sample_id") %>% 
  dplyr::mutate(rpm = count / (library_size / 10e5)) %>% 
  dplyr::filter((alignment_length >= 18) & (alignment_length <= 33)) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "WT")) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = "WT")) %>% 
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(sample_id))) %>% 
  dplyr::arrange(sample_id) 

# calculate mean and SD
read_length_stats <- 
  read_length_tb %>% 
  dplyr::group_by(genotype, alignment_length) %>% 
  dplyr::summarise(N = n(), 
                   mean = mean(rpm), 
                   median = median(rpm),
                   sd = sd(rpm), 
                   stand_err = sd / sqrt(N)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(confidence_interval = stand_err * qt((0.95 / 2) + 0.5, N - 1)) %>% 
  dplyr::arrange(desc(genotype))

# get read length distribution per sample
read_length_dist_tb <- 
  library_hist_tidy %>% 
  dplyr::select(sample_id, alignment_length, count) %>% 
  dplyr::mutate(alignment_length = as.numeric(alignment_length)) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(count_16to20nt = sum(count[alignment_length >=16 & alignment_length <= 20]),
                   count_21to23nt = sum(count[alignment_length >=21 & alignment_length <= 23]),
                   count_24to32nt = sum(count[alignment_length >=24 & alignment_length <= 32]),
                   count_33plus = sum(count[alignment_length > 32]),
                   count_total = sum(count)) %>%
  dplyr::ungroup(.) %>% 
  tidyr::pivot_longer(names_to = "length", values_to = "count", cols = -c(sample_id, count_total), names_prefix = "count_") %>% 
  dplyr::mutate(percentage = round((count / count_total), 3) * 100) %>% 
  dplyr::mutate(length = factor(length, levels = c("33plus", "24to32nt", "21to23nt", "16to20nt"), 
                                labels = c("33+", "24-32", "21-23", "16-20"))) %>% 
  dplyr::arrange(sample_id, desc(length)) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(pos = cumsum(percentage) - (0.5 * percentage)) %>%
  dplyr::ungroup(.) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "WT")) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = "WT")) %>% 
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(sample_id))) %>% 
  dplyr::arrange(sample_id)


### create plots
# # Joy Division plot
# ridge_plot <-
#   ggplot() +
#   geom_density_ridges(data = read_length_tb, 
#                       mapping = aes(x = alignment_length, y = sample_id, height = rpm, fill = sample_id), 
#                       scale = 1.5, 
#                       stat = "identity") +
#   scale_fill_viridis_d() +
#   scale_y_discrete(expand = c(0.01, 0)) +
#   scale_x_continuous(breaks = 18:33, labels = as.character(18:33)) +
#   guides(fill = FALSE) +
#   theme_ridges() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(), 
#         axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "right") 

# plot RPM
bar_plot <- 
  ggplot() +
  geom_bar(data = read_length_stats, 
           aes(x = alignment_length, y = mean, fill = genotype),
           stat = "identity",
           position = "dodge",
           width = 0.8) +
  geom_errorbar(data = read_length_stats, 
                aes(x = alignment_length, y = mean, ymin = mean - sd, ymax = mean + sd, color = genotype), 
                position = "dodge",
                size = 1, 
                width = 0.8) +
  geom_errorbar(data = read_length_stats, 
                aes(x = alignment_length, y = mean, ymin = mean, ymax = mean, color = genotype), 
                position = "dodge",
                size = 1,
                width = 0.8) +
  scale_x_discrete(limits = 18:33, breaks = 18:33, labels = as.character(18:33)) +
  scale_y_continuous(labels = fancy_scientific) + 
  scale_color_manual(values = rep("gray50", length(unique(read_length_stats$genotype)))) +
  scale_fill_viridis_d() +
  guides(color = FALSE, 
         fill = guide_legend(reverse = FALSE)) + 
  labs(subtitle = "RPM per read length (average per genotype)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.subtitle = element_text(size = 18), 
        legend.position = "bottom", 
        legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  theme(strip.text.y.right = element_text(angle = 90, size = 16))


# plot read percentage
read_length_plot <- 
  ggplot() +
  geom_bar(data = read_length_dist_tb, 
           mapping = aes(x = sample_id, y = percentage, fill = length), 
           width = 0.9, 
           stat = "identity") +
  geom_text(data = read_length_dist_tb,
            mapping = aes(x = sample_id, y = pos, label = str_c(percentage, "%")), 
            size = 6, 
            stat = "identity") +
  geom_text(data = read_length_dist_tb[(read_length_dist_tb$length == "33+" | read_length_dist_tb$length == "24-32"), ],
            mapping = aes(x = sample_id, y = pos, label = str_c(percentage, "%")), 
            size = 6, 
            stat = "identity", 
            color = "white") +
  scale_fill_viridis_d() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(subtitle = "Read length fraction per sample") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.subtitle = element_text(size = 18), 
        legend.position = "bottom", 
        legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) + 
  theme(plot.margin = unit(c(0.5, 0.5, 10, 0.5), "cm")) 


### grid and save
grid <- plot_grid(bar_plot, read_length_plot,
                  ncol = 2, 
                  align = "v", 
                  axis = "b")

# # save as pdf
# ggsave(filename = file.path(outpath, "read_width_distribution.pdf"),
#        plot = grid,
#        width = 20,
#        height = 20)

# save as object
saveRDS(grid, file.path(outpath, "grid_1.RDS"))

