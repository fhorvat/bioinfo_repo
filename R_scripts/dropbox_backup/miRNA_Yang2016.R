### INFO: 
### DATE: 18. 07. 2017.  
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
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/smallRNA/Yang_2016_PRJNA257532"
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Yang_2016_development_smallRNA_PRJNA257532/Data/Mapped/STAR_mm10"
bam_list <- list.files(path = bam_path, pattern = "s_oocyte.*bam$|s_1cell.*bam$", full.names = T)
log_list <- list.files(path = bam_path, pattern = "s_oocyte.*Log.final.out$|s_1cell.*Log.final.out$", full.names = T)

repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz"
miRNA_path <- "/common/WORK/fhorvat/reference/mouse/mm10/miRBase/miRNA_mm_20170725.gff.gz"
ensembl_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# sample table
sample_table <- 
  tibble(sample = str_replace_all(bam_list, "\\/.*\\/|_Aligned.sortedByCoord.out.bam", ""), 
         stage = str_replace_all(sample, "_r[0-9]", ""),
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
# ### read coordinates of mapped reads from .bam
# # create list of genomic coordinates 21-23 bp long
# smallRNA_list <- lapply(bam_list, function(bam_file){
# 
#   # read bam file as data frame
#   param <- Rsamtools::ScanBamParam(what = c("rname", "strand", "pos", "seq", "cigar"))
#   bam_in <- Rsamtools::scanBam(bam_file, param = param)
# 
#   # get aligned sequence
#   seq_aligned_width <-
#     GenomicAlignments::sequenceLayer(x = bam_in[[1]]$seq, cigar = bam_in[[1]]$cigar) %>%
#     width(.)
# 
#   # convert to data.frame, filter for miRNA clusters
#   bam_in[[1]]$seq <- NULL
#   bam_in[[1]]$cigar <- NULL
#   bam_gr <-
#     dplyr::bind_cols(bam_in) %>%
#     dplyr::mutate(seq_width = seq_aligned_width) %>%
#     dplyr::mutate(end = (pos + seq_width - 1),
#                   full_pos = str_c(rname, ":", pos, "-", end, "|", strand)) %>%
#     dplyr::filter(!duplicated(full_pos)) %>%
#     dplyr::select(-seq_width, seqnames = rname, start = pos, end, strand, full_pos) %>%
#     makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
#     GenomicRanges::reduce(.) %>%
#     .[width(.) >= 21 & width(.) <= 23]
# 
#   return(bam_gr)
# 
# })
# 
# # get set of miRNA ranges for all 6 samples
# bam_gr_all <- 
#   Reduce(function(gr1, gr2) c(gr1, gr2), smallRNA_list) %>% 
#   unique(.) %>% 
#   GenomicRanges::reduce(.) %>% 
#   .[width(.) >= 21 & width(.) <= 23]
# 
# # summarizeOverlaps
# bamfiles <- Rsamtools::BamFileList(bam_list, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = bam_gr_all,
#                                            reads = bamfiles,
#                                            mode = "IntersectionStrict",
#                                            singleEnd = TRUE,
#                                            ignore.strand = FALSE)
# 
# # get assay data.frame, gather to long data.frame
# assay_df <- 
#   assay(se) %>% 
#   as_tibble(.) %>% 
#   magrittr::set_colnames(., str_replace(colnames(.), "_Aligned.sortedByCoord.out.bam", "")) %>% 
#   dplyr::mutate(full_pos = str_c(seqnames(bam_gr_all), ":", start(bam_gr_all), "-", end(bam_gr_all), "|", strand(bam_gr_all))) %>% 
#   tidyr::gather(key = sample, value = count, -full_pos)


### counts and RPMs
# read assay data.frame from .RDS
assay_df <- readRDS(file = file.path(outpath, "smallRNA_reduced_assay_df.RDS"))

# get positions supported by at least 5 reads (average of replicates) in oocyte/1cell stage
smallRNA_pos <- 
  assay_df %>%
  dplyr::mutate(stage = str_replace(sample, "_r[1-9]", "")) %>% 
  dplyr::group_by(full_pos, stage) %>% 
  dplyr::summarise(count = round(mean(count), 2)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::spread(key = stage, value = count) %>% 
  dplyr::filter_at(vars(starts_with("s_")), all_vars(. >= 5)) %>% 
  dplyr::select(full_pos) %T>% 
  readr::write_csv(., path = file.path(outpath, "Yang_smallRNA_clusters_5reads.csv")) %$%
  full_pos

# get RPM values
rpm_df <- 
  assay_df %>% 
  dplyr::filter(full_pos %in% smallRNA_pos) %>%
  dplyr::left_join(., sample_table, by = "sample") %>% 
  dplyr::mutate(library_size = (library_size / 1E6), 
                rpm = (count / library_size)) %>% 
  dplyr::select(-count, -library_size, -sample) %>% 
  dplyr::group_by(full_pos, stage) %>%
  dplyr::summarise(rpm = round(mean(rpm), 3)) %>%
  dplyr::ungroup(.) %>%
  tidyr::spread(key = stage, value = rpm) %>% 
  dplyr::select(1, 3, 2)

# make GenomicRanges
smallRNA_gr <- 
  rpm_df %>% 
  tidyr::separate(full_pos, into = c("seqnames", "coordinates", "strand"), sep = ":|\\|", remove = F) %>% 
  tidyr::separate(coordinates, into = c("start", "end"), sep = "-") %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)


### overlaps
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
smallRNA_plot <- 
  smallRNA_overlaps_df %>% 
  data.table::setnames(., old = 2:3, new = c("RPM_oocyte", "RPM_1cell")) %>% 
  dplyr::filter(RPM_oocyte > 0 & RPM_1cell > 0) %>%
  dplyr::mutate(type = ifelse(!is.na(repeat_class), yes = "repeat", no = "other"), 
                type = replace(x = type, list = (!is.na(miRNA_ID) & type == "other"), values = "miRNA"), 
                type = factor(type, levels = c("other", "repeat", "miRNA"))) %>% 
  dplyr::select(1:3, 8, 4:7) %T>%
  readr::write_csv(., path = file.path(outpath, "Yang_RPM.csv")) %>% 
  dplyr::mutate(subtype = ifelse(!is.na(repeat_class), yes = repeat_class, no = "other"), 
                subtype = replace(x = subtype, list = (!is.na(miRNA_ID) & subtype == "other"), values = "miRNA"),
                subtype = factor(subtype, levels = c("other", "LINE", "LTR", "SINE", "miRNA")), 
                log2RPM_GV = log2(RPM_oocyte), 
                log2RPM_1C = log2(RPM_1cell),
                log2FC_1C_GV = (log2RPM_1C - log2RPM_GV)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(mean_RPM = mean(c(RPM_1cell, RPM_oocyte))) %>% 
  dplyr::arrange(type, subtype) %>% 
  dplyr::mutate(type = factor(type, levels = c("repeat", "miRNA", "other")), 
                subtype = factor(subtype, levels = c("SINE", "LTR", "LINE", "miRNA", "other"))) %>% 
  dplyr::filter(stringr::str_detect(string = full_pos, pattern = "chr11:115125774-115125796") | type == "miRNA") 

# MA plot
ggplot(data = smallRNA_plot, aes(x = mean_RPM, y = log2FC_1C_GV, color = type)) + 
  geom_point(size = 0.2) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000)) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_colour_manual(values = c(`repeat` = "red2", miRNA = "black", other = "gray60")) +
  xlab("mean RPM") + 
  ylab("log2FoldChange 1C/GV") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "Yang_MAplot_1CvsGV.png"), width = 9, height = 7)

# scatter plot
ggplot(data = smallRNA_plot, aes(x = log2RPM_GV, y = log2RPM_1C, color = type)) + 
  geom_point(size = 0.5) +
  scale_colour_manual(values = c(`repeat` = "red2", miRNA = "black", other = "gray60")) +
  xlab("GV log2 RPM") + 
  ylab("1C log2 RPM") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "Yang_miRNA_scatterPlot_1CvsGV.png"), width = 7, height = 7)

### color = repeat class
# MA plot
ggplot(data = smallRNA_plot, aes(x = mean_RPM, y = log2FC_1C_GV, color = subtype)) +
  geom_point(size = 0.2) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000)) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_colour_manual(values = c(SINE = "red2", LTR = "green", LINE = "blue", miRNA = "black", other = "gray60")) +
  xlab("mean RPM") +
  ylab("log2FoldChange 1C/GV") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "Yang_MAplot_color_1CvsGV.png"), width = 9, height = 7)

# scatter plot
ggplot(data = smallRNA_plot, aes(x = log2RPM_GV, y = log2RPM_1C, color = subtype)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = c(SINE = "red2", LTR = "green", LINE = "blue", miRNA = "black", other = "gray60")) +
  xlab("GV log2 RPM") +
  ylab("1C log2 RPM") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "Yang_scatterPlot_color_1CvsGV.png"), width = 7, height = 7)