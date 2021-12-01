library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files")

################################################################################## functions
# Clustal Omega (vfranke)
ClustalO <- function(infile, outfile, ClustalO = NULL, threads = 12, 
                     what = 'DNA', force = T, format = 'fa'){
  
  if(!what %in% c('AA', 'DNA'))
    stop('can only align DNA or AA')
  if(is.null(ClustalO))
    ClustalO <- '/common/WORK/fhorvat/programi/clustal-omega-1.2.3/bin/clustalo'
  
  ### checks whether the file exists and whether to force the outfile
  ### if the file does exist and the force is off he reads the file
  if((file.exists(outfile) & force == T) | !file.exists(outfile)){
    cat('Running the alignmnent...\n')
    threads <- paste('--threads=', threads, sep = '')
    format <- paste('--outfmt=', format, sep = '')
    command <- paste(ClustalO, '-i', infile, '-o', outfile,'--force' ,threads, format)
    system(command)
  }
  
  cat('Returning the results...\n')
  if(what == 'AA')
    a <- readAAStringSet(outfile, format = 'fasta')
  if(what == 'DNA')
    a <- readDNAStringSet(outfile, format = 'fasta')
  
  return(a)
}

# get consensus sequence from repbase
getLTRConsensus <- function(LTR_name){
  
  # get lines from repbase library between LTR name and next LTR element, write as .txt
  command_consensus <- paste0("cat /common/WORK/fhorvat/programi/repeatMasker/Repbase_library/Libraries/RepeatMaskerLib.embl", 
                              " | ", "sed -n '/", LTR_name, "/,/^ID/p'", 
                              " > ", 
                              "/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, ".txt")
  system(command_consensus)
  
  # get sequence consensus
  ltr_consensus <- 
    read_lines(paste0(LTR_name, ".txt")) %>% 
    .[(grep("SQ", .) + 1) : (grep("//", .) - 1)] %>% 
    str_replace_all("[0-9]*|", "") %>% 
    str_replace_all(" ", "") %>%
    paste(collapse = "") %>%
    toupper()
  
  # removes .txt wit LTR data
  system(paste0("rm ", "/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, ".txt")) 
  
  return(ltr_consensus)
}

countAndFPKM <- function(features_granges, singleEndBoolean, sample_df){
  
  # counts
  register(MulticoreParam())
  se <- summarizeOverlaps(features = features_granges, 
                          reads = BamFileList(sample_df$track_path, yieldSize = 2000000), 
                          mode = "Union", 
                          singleEnd = singleEndBoolean, 
                          ignore.strand = TRUE)
  
  # FPKM
  fpkm_df <- 
    as.data.frame(assay(se)) %>% 
    set_colnames(sample_df$sample_name) %>% 
    mutate(width = width(features_granges), 
           fullName = features_granges$fullName)
  
  invisible(lapply(X = sample_df$sample_name,
                   FUN = function(X) fpkm_df[, X] <<- fpkm_df[, X] / (sample_df[sample_df$sample_name == X, "library_size"] * (fpkm_df$width / 1000))))
  
  return(fpkm_df)
}

alignLTRClustalOAndPlot <- function(rptmsk_LTR, LTR_set_name){
  
  # get sequences 
  LTR_seq <- as.character(getSeq(x = Mmusculus, rptmsk_LTR))
  
  # merge consensus sequence with the consensus
  LTR_seq <- c(ltr_consensus, LTR_seq)
  
  # write sequences as FASTA
  write.fasta(as.list(LTR_seq),
              nbchar = 80,
              names = c("consensus", paste0(rptmsk_LTR$fullName)),
              as.string = TRUE,
              file.out = paste0(LTR_name, "_", LTR_set_name, ".fasta"),
              open = "w")
  
  # get sequences aligned by Clustal Omega
  LTR_seq_aligned <- ClustalO(infile = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, ".fasta"), 
                              outfile = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, "_aligned.fasta"), 
                              threads = 8, 
                              force = T)
  
  # melt data.frame for plot
  LTR_seq_aligned_df <- 
    as.data.frame(as.matrix(LTR_seq_aligned), stringsAsFactors = F) %>% 
    mutate(ID = rownames(.))
    
  if(LTR_set_name == "expressed"){
    LTR_seq_aligned_df %<>%
      mutate(ID = factor(ID, rev(c("consensus", fpkm_LTR$fullName))))
  }
  
  LTR_seq_aligned_df %<>% 
    dplyr::select(extract(., 1, ) %>% as.character() %>% equals("-") %>% not() %>% which()) %>% 
    melt(id.vars = "ID") %>% 
    mutate(value = replace(value, value == "-", NA), 
           value = factor(value, levels = c("A", "C", "G", "T")))
  
  # plot
  ggplot(LTR_seq_aligned_df, aes(x = variable, y = ID)) +
    geom_tile(aes(fill = value)) +
    scale_fill_manual(values = color_pallete, 
                      breaks = c("A", "C", "G", "T")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(),
          # axis.text.y = element_text(size = 1, hjust = 1, vjust = 0),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggsave(paste0(LTR_name, "_", LTR_set_name, "_msa.pdf"))
  
  return(paste0(LTR_name, "_", LTR_set_name))
  
}

################################################################################## reading data
# repeatMasker
rptmsk <- read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t")

# Vedran's figure 3 data 
repeats_list <- readRDS("/common/WORK/vfranke/Projects/PSvoboda_MT/Results/MT_FindRepeats/mm/160907.mm.Results.rds", refhook = NULL)

# mannualy set colors to the stages
# (A) green; (C) blue; (G) yellow; (T) red
color_pallete <- c("#66CD00", "#104E8B", "#FFD700", "#CD2626")
names(color_pallete) <- c("A", "C", "G", "T")

# Fugaku data samples
sample_df_LTR <- 
  data.frame(track_path = list.files(path = "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2", 
                                     pattern = "*.bam$", 
                                     recursive = T, 
                                     full.names = T), 
             log_path = list.files(path = "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2", 
                                   pattern = "*Log.final.out", 
                                   recursive = T, 
                                   full.names = T),
             stringsAsFactors = F) %>%
  mutate(sample_name = gsub("^/.*/|\\.bam", "", track_path), 
         library_size = sapply(X = log_path, function(X) as.integer(read.delim(X, header = F, stringsAsFactors = F)[8, 2]) / 10e6)) %>% 
  dplyr::select(sample_name, track_path, library_size) %>%
  right_join(data.frame(sample_name = "s_GV.WE", stringsAsFactors = F), by = "sample_name")

################################################################################### get coordinates of specific LTR element (expressed and random)
# set LTR name
LTR_name <- "MTA_Mm"

# get consensus sequence from repbase
ltr_consensus <- getLTRConsensus(LTR_name)

# full contribution to a 5’ exon in protein-coding genes
selected_MT <- 
  repeats_list$RepeatsSelected_MALR %>% 
  filter(ex.category == "5' exon", 
         ex.category.complete == "Complete", 
         gene_biotype == "protein_coding", 
         repName == LTR_name) %>% 
  mutate(fullName = paste0(rep.coord, ":", 
                           rep.strand, "|", 
                           repName))

# find expressed LTR in repeatMasker
rptmsk_LTR_expressed <-
  rptmsk %>% 
  filter(grepl(LTR_name, element_name) & !grepl("int", element_name)) %>%
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name)) %>% 
  dplyr::select(-element_class) %>% 
  filter(fullName %in% selected_MT$fullName) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# get random LTRs (the same number as above), take only those with length = consensus length +- 5%
set.seed(1111)
rptmsk_LTR_random <-
  rptmsk %>% 
  filter(grepl(LTR_name, element_name) & !grepl("int", element_name)) %>%
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name), 
         width = end - start + 1) %>% 
  filter(width > nchar(ltr_consensus) - (0.05 * nchar(ltr_consensus)), 
         width < nchar(ltr_consensus) + (0.05 * nchar(ltr_consensus))) %>% 
  dplyr::select(-element_class, -width) %>% 
  sample_n(nrow(selected_MT)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

################################################################################### calculate mean FPKM expression in 3' half of LTR in MII/GV, order 
# 3' half of LTR
rptmsk_LTR_expressed_3prime_half <- GenomicRanges::resize(rptmsk_LTR_expressed, width = round(width(rptmsk_LTR_expressed) / 2), fix = "end")

# FPKM ordered by expression in GV
fpkm_LTR <- countAndFPKM(features = rptmsk_LTR_expressed_3prime_half, singleEndBoolean = F, sample_df = sample_df_LTR)  
write_csv(fpkm_LTR, paste0(LTR_name, "_expressed_FPKM_in_3PrimeHalf_GV_test.csv"))

################################################################################### align with ClustalO and plot
alignLTRClustalOAndPlot(rptmsk_LTR = rptmsk_LTR_expressed, LTR_set_name = "expressed")
alignLTRClustalOAndPlot(rptmsk_LTR = rptmsk_LTR_random, LTR_set_name = "random")

