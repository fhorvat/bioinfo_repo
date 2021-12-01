### INFO: 
### DATE: 10. 08. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()[!(ls() %in% c("rptmsk_retro_gr", "miRNA_gr"))]); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/smallRNA/Yang_2016_PRJNA257532")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/smallRNA"
table_list <- list.files(path = outpath, pattern = "*smallRNA_clusters_5reads.csv", full.names = T, recursive = T)

bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq"
bam_list <- 
  list.files(path = bam_path, pattern = "s_.*bam$", full.names = T, recursive = T) %>% 
  stringr::str_subset(., pattern = "oocyte|1cell|zygote")
log_list <- 
  bam_list %>% 
  stringr::str_replace(string = ., pattern = "_Aligned.sortedByCoord.out.bam", replacement = "_Log.final.out") 

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz"
ensembl_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# sample table
sample_table <- 
  tibble(sample = str_c(str_replace_all(bam_list, "\\/.*\\/|_Aligned.sortedByCoord.out.bam", "")),
         data = str_extract(bam_list, "Yang|GarciaLopez"), 
         log_path = log_list) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(library_size = readr::read_lines(file = log_path) %>% 
                  magrittr::extract(9) %>% 
                  str_extract(., "[0-9].*") %>% 
                  as.integer(.)) %>% 
  dplyr::select(-log_path)

# repeatMaskerVIZ - retrotransposons only
rptmsk_retro_gr <-
  read_delim(file = repeatmasker_path, delim = "\t") %>%
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>%
  dplyr::filter(repClass %in% c("LINE", "SINE", "LTR")) %>%
  dplyr::mutate(full_pos = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repClass, "|", repName)) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# miRNA from ENSEMBL
miRNA_gr <-
  read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>%
  GffToGRanges(., filter = "exon") %>%
  as.data.frame(.) %>%
  as_tibble(.) %>%
  dplyr::filter(gene_biotype == "miRNA") %>%
  dplyr::select(seqnames, start, end, strand, gene_id) %>%
  dplyr::mutate(full_pos = str_c(seqnames, ":", start, "-", end, "|", strand, "|", gene_id)) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## MAIN CODE
### overlaping position by ranges
# get ranges
GL_ranges <- 
  readr::read_csv(file = table_list[str_detect(table_list, "GarciaLopez")]) %>% 
  tidyr::separate(full_pos, into = c("seqnames", "coordinates", "strand"), sep = ":|\\|") %>% 
  tidyr::separate(coordinates, into = c("start", "end"), sep = "-") %>% 
  GenomicRanges::makeGRangesFromDataFrame(.)

Yang_ranges <- 
  readr::read_csv(file = table_list[str_detect(table_list, "Yang")]) %>% 
  tidyr::separate(full_pos, into = c("seqnames", "coordinates", "strand"), sep = ":|\\|") %>% 
  tidyr::separate(coordinates, into = c("start", "end"), sep = "-") %>% 
  GenomicRanges::makeGRangesFromDataFrame(.)

# create union of both ranges
ranges_all <- 
  c(GL_ranges, Yang_ranges) %>% 
  unique(.) %>%
  GenomicRanges::reduce(.) %>%
  .[width(.) >= 21 & width(.) <= 23]

### counts and RPMs
# # summarizeOverlaps
# bamfiles <- Rsamtools::BamFileList(bam_list, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = ranges_all,
#                                            reads = bamfiles,
#                                            mode = "IntersectionStrict",
#                                            singleEnd = TRUE,
#                                            ignore.strand = FALSE)

se <- readRDS(file = file.path(outpath, "sumOverlaps_Yang_GL.RDS"))

# get assay data.frame, gather to long data.frame
assay_df <-
  assay(se) %>%
  as_tibble(.) %>%
  magrittr::set_colnames(., str_replace(colnames(.), "_Aligned.sortedByCoord.out.bam", "")) %>%
  dplyr::mutate(full_pos = str_c(seqnames(ranges_all), ":", start(ranges_all), "-", end(ranges_all), "|", strand(ranges_all))) %>%
  tidyr::gather(key = sample, value = count, -full_pos) 

# get RPM values
rpm_df <- 
  assay_df %>% 
  dplyr::left_join(., sample_table, by = "sample") %>% 
  dplyr::mutate(library_size = (library_size / 1E6), 
                rpm = round((count / library_size), 3)) %>% 
  dplyr::mutate(sample = stringr::str_c(data, "_", sample) %>% 
                  stringr::str_replace_all(string = ., pattern = "s_|_r[1-9]", replacement = "")) %>% 
  dplyr::select(-count, -library_size, -data) %>% 
  dplyr::group_by(full_pos, sample) %>% 
  dplyr::summarise(rpm = round(mean(rpm), 3)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::spread(key = sample, value = rpm) 
  

### overlaps
# make GenomicRanges
smallRNA_gr <- 
  rpm_df %>% 
  tidyr::separate(full_pos, into = c("seqnames", "coordinates", "strand"), sep = ":|\\|", remove = F) %>% 
  tidyr::separate(coordinates, into = c("start", "end"), sep = "-") %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# overlap with repeatMasker
smallRNA_retro_overlaps <- GenomicAlignments::findOverlaps(smallRNA_gr, rptmsk_retro_gr, ignore.strand = F, minoverlap = min(width(smallRNA_gr)))
smallRNA_retro_df <- tibble(full_pos = smallRNA_gr$full_pos[queryHits(smallRNA_retro_overlaps)],
                            repeat_class = rptmsk_retro_gr$repClass[subjectHits(smallRNA_retro_overlaps)], 
                            repeat_family = rptmsk_retro_gr$repFamily[subjectHits(smallRNA_retro_overlaps)], 
                            repeat_name = rptmsk_retro_gr$repName[subjectHits(smallRNA_retro_overlaps)])

# overlap with miRNA
smallRNA_miRNA_overlaps <- GenomicAlignments::findOverlaps(smallRNA_gr, miRNA_gr, ignore.strand = F, minoverlap = min(width(smallRNA_gr)))
smallRNA_miRNA_df <- tibble(full_pos = smallRNA_gr$full_pos[queryHits(smallRNA_miRNA_overlaps)], 
                            miRNA_ID = miRNA_gr$gene_id[subjectHits(smallRNA_miRNA_overlaps)])

# join small RNA reads with repeatMasker and miRNA overlaps
smallRNA_overlaps_df <- dplyr::left_join(rpm_df, full_join(smallRNA_retro_df, smallRNA_miRNA_df, by = "full_pos"), by = "full_pos") 


### plot, write table
# plot table
rpm_plot <- 
  smallRNA_overlaps_df %>% 
  dplyr::filter_at(vars(matches("Yang.*|GarciaLopez.*")), all_vars(. > 0)) %>%
  dplyr::mutate_at(vars(matches("Yang.*|GarciaLopez.*")), all_vars(log2(.))) %>% 
  dplyr::mutate(type = ifelse(!is.na(repeat_class), yes = "repeat", no = "other"), 
                type = replace(x = type, list = (!is.na(miRNA_ID) & type == "other"), values = "miRNA"), 
                type = factor(type, levels = c("other", "repeat", "miRNA")),
                subtype = ifelse(!is.na(repeat_class), yes = repeat_class, no = "other"), 
                subtype = replace(x = subtype, list = (!is.na(miRNA_ID) & subtype == "other"), values = "miRNA"),
                subtype = factor(subtype, levels = c("other", "LINE", "LTR", "SINE", "miRNA")), 
                log2FC_1C_GV_GarciaLopez = (GarciaLopez_zygote - GarciaLopez_oocyte), 
                log2FC_1C_GV_Yang = (Yang_1cell - Yang_oocyte)) %>% 
  dplyr::arrange(type, subtype) %>% 
  dplyr::mutate(type = factor(type, levels = c("repeat", "miRNA", "other")), 
                subtype = factor(subtype, levels = c("SINE", "LTR", "LINE", "miRNA", "other")))

# plot one color
ggplot(data = rpm_plot, aes(x = log2FC_1C_GV_Yang, y = log2FC_1C_GV_GarciaLopez, color = type)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_colour_manual(values = c(`repeat` = "red2", miRNA = "black", other = "gray60")) +
  xlab("log2FC 1C/GV Yang") + 
  ylab("log2FC 1C/GV GarciaLopez") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "log2FC_1C_GV_GarciaLopez_vs_Yang.png"), width = 7, height = 7)


# plot 3 colors
ggplot(data = rpm_plot, aes(x = log2FC_1C_GV_Yang, y = log2FC_1C_GV_GarciaLopez, color = subtype)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_colour_manual(values = c(SINE = "red2", LTR = "green", LINE = "blue", miRNA = "black", other = "gray60")) +
  xlab("log2FC 1C/GV Yang") + 
  ylab("log2FC 1C/GV GarciaLopez") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "log2FC_1C_GV_GarciaLopez_vs_Yang_color.png"), width = 7, height = 7)

