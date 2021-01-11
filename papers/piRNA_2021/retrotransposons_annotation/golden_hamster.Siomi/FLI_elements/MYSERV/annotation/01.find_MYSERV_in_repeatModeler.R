### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV/annotation")

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
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)

library(BSgenome.Maur.UCSC.Siomi)
library(seqinr)
library(Biostrings)
library(systemPipeR)
library(stringdist)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

### repeatModeler
# joined repeatModeler path
rmod_path <- file.path(genome_dir, "RepeatModeler/RepeatMasker", "rmsk.Siomi.20200728.RepeatModeler.joined_rmsk_id.fa.out.gz")

# clean repeatModeler path
rmod_clean_path <- file.path(genome_dir, "RepeatModeler/RepeatMasker", "rmsk.Siomi.20200728.RepeatModeler.clean.fa.out.gz")

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

# read joined repeatModeler
rmod_tb <- readr::read_delim(rmod_path, delim = "\t")

# read clean repeatModeler
rmod_clean <- readr::read_delim(rmod_clean_path, delim = "\t")

######################################################## MAIN CODE
### clean data
# get widths
rmsk_tb %<>%    
  # dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "+", strand)) %>%
  # GRanges(.) %>%
  # as_tibble(.) %>%
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

# get repeatModeler LTRs as table
rmod_ltr_tb <- 
  rmod_clean %>% 
  dplyr::filter(repClass == "LTR")

# get repeatModeler LTRs as GRanges
rmod_ltr_gr <-
  rmod_ltr_tb %>% 
  GRanges(.)


### find repeatModeler annotation
# filter MYSERV6-int
rmsk_tb_filt <-
  rmsk_clean %>%
  dplyr::filter(repName %in% c("MYSERV6-int", "MYSERV-int", "MYSERV16_I")) %>% 
  GRanges(.)

# overlap with repeatModeler
overlap <- findOverlaps(rmsk_tb_filt, rmod_ltr_gr, ignore.strand = F)

# get repeatModeler hits
rmod_filt <- 
  rmod_ltr_gr[subjectHits(overlap)] %>% 
  as_tibble(.) %>% 
  distinct(.) %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::arrange(-count)

# get joined repeatModeler
rmod_joined_filt <- 
  rmod_tb %>% 
  dplyr::filter(repName %in% c("ltr-1_family-5", "ltr-1_family-6", "ltr-1_family-16")) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1) %>% 
  dplyr::arrange(-width) %>% 
  dplyr::filter(insertion_class != "interupted")

# plot histogram of widths
hist_width <-
  ggplot(data = rmod_joined_filt, aes(width, fill = repName)) +
  geom_histogram(binwidth = 10) +
  scale_x_continuous(limits = c(0, 10000), breaks = seq(0, 10000, 1000)) +
  facet_grid(rows = vars(repName), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, str_c("MYSERV", "repeatModeler", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 5)


### plot only LTR family 16
# get joined repeatModeler
ltr_family_16_tb <-
  rmod_joined_filt %>%
  dplyr::filter(repName == c("ltr-1_family-16"))

# plot histogram of widths
hist_width <-
  ggplot(data = ltr_family_16_tb, aes(width, fill = repName)) +
  geom_histogram(binwidth = 5) +
  scale_x_continuous(limits = c(2000, 3000), breaks = seq(2000, 3000, 100)) +
  # facet_grid(rows = vars(repName), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, str_c("ltr-1_family-16", "repeatModeler", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 5)


### filter based on histogram of widths
# filter everything between 2700-3000
ltr_family_16_tb %<>% 
  dplyr::filter(width >= 2700, width <= 3000)

### save table and fasta
# save 
readr::write_csv(ltr_family_16_tb, file = file.path(outpath, str_c("MYSERV.ltr-1_family-16.2.7_to_3k", "csv", sep = ".")))

# get DNAStringSet
ltr_family_16_seq <- 
  ltr_family_16_tb %>% 
  GRanges(.) %>% 
  getSeq(BSgenome.Maur.UCSC.Siomi, .)

# set names
names(ltr_family_16_seq) <- str_c("MYSERV", ltr_family_16_tb$rmsk_id, sep = ".")

# save as fasta
Biostrings::writeXStringSet(x = ltr_family_16_seq, 
                            filepath = file.path(outpath, str_c("MYSERV.ltr-1_family-16.2.7_to_3k", "fasta", sep = ".")))




