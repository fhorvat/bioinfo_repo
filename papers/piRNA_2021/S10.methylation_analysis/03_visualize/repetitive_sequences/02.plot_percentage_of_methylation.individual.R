### INFO: 
### DATE: Wed Dec 09 18:53:14 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/methylation/repetitive")

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

library(GenomicRanges)

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

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# get features path
features_path <- file.path(inpath, "rmsk.chosen_LTRs_and_L1s.GRanges.RDS")

# get methylation calls path
meth_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi.trimmed"
meth_path <- list.files(meth_path, ".*\\.PE_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz", full.names = T)

######################################################## READ DATAl
# read table coverage data
meth_tb <- purrr::map(meth_path, function(path){
  
  # read table
  readr::read_delim(path, delim = "\t", col_names = c("seqnames", "start", "end", "percentage_meth", "count_meth", "count_unmeth")) %>% 
    dplyr::mutate(sample_id = basename(path) %>% str_remove(., "_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz|\\.merged\\.bismark\\.cov\\.gz"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO"))) %>% 
  tidyr::unite(col = coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(coordinates, count_meth, count_unmeth, genotype)

# read features
features_list <- readRDS(features_path)
features_list <- features_list[!names(features_list) %in% c("SINEB2", "RLTR3x ERV1", "L1Mus")]

######################################################## MAIN CODE
# set minimal count
count_min <- 4

# create hierarchy table: Lx5, L1Mus, MuLV, MYSERVx, RLTR3x ERV1, IAP & ORR1
hierch_tb <- tibble(hit_class = c("Lx5", "MuLV", "MYSERVx", "IAP", "ORR1"),
                    hierch_number = 1:5)

# calculate total number of reads per site
meth_sites <- 
  meth_tb %>% 
  dplyr::mutate(count = count_meth + count_unmeth)
  

### get methylation levels for different genotypes separately
meth_sites_genotype <- purrr::map(c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO"), function(one_genotype){
  
  # filter table, 
  meth_sites_tb <- 
    meth_sites %>% 
    dplyr::filter(genotype == one_genotype, 
                  count >= count_min) %>% 
    tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ", remove = F)
  
  # get methylation sites as GRanges
  meth_sites_gr <- GRanges(meth_sites_tb)
  
  # overlap with features
  overlaps <- findOverlaps(meth_sites_gr, features_list, ignore.strand = T)
  
  # create table with hits
  meth_sites_hits <- 
    tibble(hit_coordinate = meth_sites_gr[queryHits(overlaps)]$coordinates,
           hit_class = names(features_list)[subjectHits(overlaps)]) %>% 
    dplyr::left_join(., hierch_tb, by = "hit_class") %>% 
    dplyr::group_by(hit_coordinate) %>% 
    dplyr::filter(hierch_number == min(hierch_number)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::select(hit_coordinate, hit_class) %>% 
    dplyr::left_join(., meth_sites_tb, by = c("hit_coordinate" = "coordinates")) %>% 
    dplyr::mutate(meth_percentage = 100*(count_meth / (count_meth + count_unmeth)))
  
  # return
  return(meth_sites_hits)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::select(hit_coordinate, hit_class, meth_percentage, count_meth, count_unmeth, genotype)


### plot barplot
# plot table
plot_tb <- 
  meth_sites_genotype %>% 
  dplyr::group_by(hit_class, genotype) %>% 
  dplyr::summarise(count_meth = sum(count_meth), 
                   count_unmeth = sum(count_unmeth), 
                   meth_percentage = 100*(count_meth / (count_meth + count_unmeth))) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(hit_class = factor(hit_class, levels = hierch_tb$hit_class))

# save as table
readr::write_csv(plot_tb, file.path(outpath, str_c("methylation_percentage", "repeat_classes", "individual_coverage_filtering", "csv", sep = ".")))

# visualize
bar_plot <-
  ggplot(plot_tb, aes(x = hit_class, y = meth_percentage, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.8) +
  scale_fill_manual(values = c("Mov10l1_WT" = "gray80",  "Mov10l1_HET" = "gray60", "Mov10l1_KO" = "black")) +
  ylab("percentage of methylation") +
  xlab("repeat class") +
  # guides(fill = TRUE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# save
ggsave(plot = bar_plot,
       filename = file.path(outpath, str_c("barplot", "methylation_percentage", "repeat_classes", "individual_coverage_filtering", "png", sep = ".")),
       width = 12, height = 10)



# ### plot violin plots
# # plot table
# plot_tb_violin <- 
#   meth_sites_genotype %>% 
#   dplyr::ungroup(.) %>% 
#   dplyr::mutate(hit_class = factor(hit_class, levels = hierch_tb$hit_class))
# 
# # visualize
# violin_plot <-
#   ggplot(plot_tb_violin, aes(x = hit_class, y = meth_percentage, fill = genotype)) +
#   geom_violin(position = "dodge", width = 0.75) +
#   geom_jitter(shape = 16, position = ggplot2::position_jitterdodge(0.1), color = "grey70") +
#   scale_fill_manual(values = c("Mov10l1_WT" = "gray80",  "Mov10l1_HET" = "gray60", "Mov10l1_KO" = "black")) +
#   ylab("percentage of methylation") +
#   xlab("repeat class") +
#   # guides(fill = TRUE) +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#         axis.text.x = element_text(size = 12, angle = 90),
#         axis.text.y = element_text(size = 10, vjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom")
# 
# # save
# ggsave(plot = violin_plot,
#        filename = file.path(outpath, str_c("violin", "methylation_percentage", "repeat_classes", "individual_coverage_filtering", "png", sep = ".")),
#        width = 12, height = 10)




