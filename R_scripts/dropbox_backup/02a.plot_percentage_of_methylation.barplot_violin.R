### INFO: 
### DATE: Wed Dec 09 18:53:14 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/methylation/non_repetitive")

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
features_path <- file.path(inpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.filt.exons_introns_promoters.RDS")

# get methylation calls path
meth_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi"
meth_path <- c(str_c(meth_path, "/s_GV_Mov10l1_KO_So821_r1.PE_bismark_bt2_pe.deduplicated.bismark.cov.gz"),
               str_c(meth_path, "/merge_HET/s_GV_Mov10l1_HET_r1.PE.merged.bismark.cov.gz"))

# meth_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi"
# meth_path <- c(str_c(meth_path, "/s_GV_Mov10l1_KO_So821_r1.PE_bismark_bt2_pe.deduplicated.bismark.cov.gz"), 
#                str_c(meth_path, "/s_GV_Mov10l1_HET_So811_r1.PE_bismark_bt2_pe.deduplicated.bismark.cov.gz"))

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

######################################################## MAIN CODE
# get sites with enough reads (4+ in both KO and HET)
meth_sites <- 
  meth_tb %>% 
  dplyr::mutate(count = count_meth + count_unmeth) %>% 
  tidyr::pivot_wider(id_cols = coordinates, values_from = "count", names_from = "genotype") %>% 
  dplyr::mutate_all(~(replace(., is.na(.), 0))) %>% 
  dplyr::filter(Mov10l1_HET >= 4 & Mov10l1_KO >= 4)

# create hierarchy table
hierch_tb <- tibble(hit_class = c("LINE", "LTR", "SINE", "rmsk_other", 
                                  "promoters_expressed", "promoters_non_expressed", 
                                  "exons", 
                                  "introns"),
                    hierch_number = 1:8)

# get methylation sites as GRanges
meth_sites_gr <- 
  meth_sites %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
  GRanges(.)

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
  dplyr::select(hit_coordinate, hit_class)

# filter methylation table, add class
meth_tb_filt <- 
  meth_tb %>% 
  dplyr::filter(coordinates %in% meth_sites$coordinates) %>% 
  dplyr::mutate(meth_percentage = 100*(count_meth / (count_meth + count_unmeth))) %>% 
  tidyr::pivot_wider(coordinates, names_from = genotype, values_from = meth_percentage) %>% 
  dplyr::left_join(., meth_sites_hits, by = c("coordinates" = "hit_coordinate")) %>% 
  dplyr::mutate(hit_class = replace(hit_class, is.na(hit_class), "intergenic"))

### test siginficance between KO/HET values
# create table with counts
test_tb <- 
  meth_tb %>% 
  dplyr::filter(coordinates %in% meth_sites$coordinates) %>%
  dplyr::left_join(., meth_sites_hits, by = c("coordinates" = "hit_coordinate")) %>% 
  dplyr::mutate(hit_class = replace(hit_class, is.na(hit_class), "intergenic")) %>% 
  dplyr::group_by(hit_class, genotype) %>% 
  dplyr::summarise(count_meth = sum(count_meth), 
                   count_unmeth = sum(count_unmeth)) %>%
  dplyr::ungroup(.) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Mov10l1_KO", "Mov10l1_HET"))) %>% 
  dplyr::arrange(hit_class, genotype)

# test for some hit class 
test_tb %>% 
  dplyr::filter(hit_class == "intergenic") %>%
  dplyr::select(-hit_class) %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "genotype") %>% 
  chisq.test(x = ., simulate.p.value = T)


### plot violin plots
# plot table
plot_tb <- 
  meth_tb_filt %>% 
  # dplyr::select(-coordinates) %>% 
  tidyr::pivot_longer(cols = -c(coordinates, hit_class), names_to = "genotype", values_to = "percentage") %>% 
  dplyr::ungroup(.) %>% 
  # dplyr::filter(hit_class != "rmsk") %>% 
  dplyr::mutate(hit_class = factor(hit_class, levels = c("promoters_expressed", "promoters_non_expressed", 
                                                         "exons", "introns", "intergenic", 
                                                         "LINE", "LTR", "SINE", "rmsk_other")))





# plot repeats
plot_list <- list(repeats = plot_tb %>% dplyr::filter(hit_class %in% c("LINE", "LTR", "SINE", "rmsk_other")), 
                  non_repetitive = plot_tb %>% dplyr::filter(!(hit_class %in% c("LINE", "LTR", "SINE", "rmsk_other"))))

# plot in loot
purrr::map(names(plot_list), function(name){
  
  # visualize
  violin_plot <-
    ggplot(plot_list[[name]], aes(x = hit_class, y = percentage, fill = genotype)) +
    geom_violin(position = "dodge", width = 0.75) +
    # geom_jitter(shape = 16, position = ggplot2::position_jitterdodge(0.1), color = "grey70") +
    scale_fill_manual(values = c("Mov10l1_HET" = "gray60", "Mov10l1_KO" = "black")) +
    ylab("percentage of methylation") +
    xlab("genomic region") +
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
  ggsave(plot = violin_plot,
         filename = file.path(outpath, str_c("violin", "methylation_percentage", name, "1kb_prom", "png", sep = ".")),
         width = 12, height = 10)
  
})


### plot barplot
# plot table
plot_tb <- 
  meth_tb_filt %>% 
  dplyr::select(-coordinates) %>% 
  tidyr::pivot_longer(cols = -hit_class, names_to = "genotype", values_to = "percentage") %>% 
  dplyr::group_by(hit_class, genotype) %>% 
  dplyr::summarise(mean_meth = mean(percentage), 
                   sd = sd(percentage), 
                   se = standard.error(percentage)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(hit_class = factor(hit_class, levels = c("promoters_expressed", "promoters_non_expressed", 
                                                         "exons", "introns", "intergenic", 
                                                         "LINE", "LTR", "SINE", "rmsk_other")))

# save as table
readr::write_csv(plot_tb, file.path(outpath, "mean_methylation_percentage_per_category.1kb_promoters.ooctye.csv"))

# split to list
plot_list <- list(repeats = plot_tb %>% dplyr::filter(hit_class %in% c("LINE", "LTR", "SINE", "rmsk_other")), 
                  non_repetitive = plot_tb %>% dplyr::filter(!(hit_class %in% c("LINE", "LTR", "SINE", "rmsk_other"))))

# plot in loot
purrr::map(names(plot_list), function(name){
  
  # visualize
  bar_plot <-
    ggplot(plot_list[[name]], aes(x = hit_class, y = mean_meth, fill = genotype)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.8) +
    geom_errorbar(aes(ymin = mean_meth - se, ymax = mean_meth + se, width = 0.8),
                  stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("Mov10l1_HET" = "gray60", "Mov10l1_KO" = "black")) +
    ylab("percentage of methylation") +
    xlab("genomic region") +
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
         filename = file.path(outpath, str_c("barplot", "methylation_percentage", name, "1kb_prom", "png", sep = ".")),
         width = 12, height = 10)
  
})

