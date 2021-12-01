### INFO: reproduces Rosa's Figure 4A from Genome Reasearch 2017 paper
### DATE: Mon Mar 11 15:47:01 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/Rosa")

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
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# bam paths
# bam_path <- "%BAM_PATH"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/Rosa/datasets/Gan_2013_NatCommun_GSE35005/filtered_reads/s_elongative_ST_r1.SE.bam"
bam_name <- basename(bam_path) %>% str_remove(., "\\.bam")

# range file path
range_path <- "/common/WORK/rosa/MT_elements/Data/MT_MT2_ORR_MLT_B1_B2_L1_IAP_elements.nogene.txt"

######################################################## READ DATA
# read in info for MT elements
mt <- read.delim(range_path, stringsAsFactors = FALSE)

######################################################## MAIN CODE
# make ranges for MT elements
mtRanges <- GRanges(IRanges(mt$genoStart, mt$genoEnd), 
                    seqnames = Rle(mt$genoName), 
                    strand = mt$strand, 
                    name = mt$repName)

# check if the experiment is paired end or single end
paired <- str_detect(bam_name, "\\.PE")

# read bam file
if (paired){
  reads <- readGAlignmentPairs(bam_path, use.names = TRUE)
} else {
  reads <- readGAlignments(bam_path, use.names = T)
}

# find all overlaps between reads and ranges
my.readOverlaps <- findOverlaps(mtRanges, reads, ignore.strand = T)

# create table with overlaps
my.readsTable <- data.frame(ranges = as.character(mtRanges$name[queryHits(my.readOverlaps)]), 
                            read_name = names(reads)[subjectHits(my.readOverlaps)])

# # remove bam file from memory
# rm(reads); gc()

# get unique reads/features 
my.uniqueReadsTable <- unique(my.readsTable)

# count reads
counts <- my.uniqueReadsTable[!(my.uniqueReadsTable$read_name %in% names(which(table(my.uniqueReadsTable$read_name) > 1))), ]

# get the final counts per repeatClass
final.counts <- as.data.frame(table(counts$ranges))

# clean and save table as .RDS
rownames(final.counts) <- final.counts[, 1]
column <- gsub("\\.bam", "", gsub("^.*\\/", "", bam_path))
final.counts[, 1] <- NULL
colnames(final.counts) <- column
save(final.counts, file = file.path(outpath, str_c(bam_name, ".read_counts.RDS")))

