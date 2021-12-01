#!/common/WORK/fhorvat/programi/R/R-3.4.3/bin/Rscript
### INFO: qsub -q MASTER -l select=ncpus=10:mem=800gb -N pbs.02.velocyto_saveRDS.pbs.R -j oe 02.velocyto_saveRDS.pbs.R
### DATE: 27. 11. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Analysis/velocyto/oocyte_no_2C")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(velocyto.R)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS
# modified read.smartseq2.bams (in parallel::mclapply mc.preschedule = F)
read.smartseq2.bams.modified <- function(bam.files, annotation.file, min.exon.count = 100, n.cores = defaultNCores()){

  cat("reading gene annotation ... ")
  x <- read.delim(annotation.file, header = F, sep = "\t", stringsAsFactors = F)
  genes <- data.frame(name = x[, 1], chr = x[, 3], strand = x[, 4], start = x[, 5], end = x[, 6], stringsAsFactors = F)
  genes$p5 <- genes$start
  genes$p5[genes$strand == "-"] <- genes$end[genes$strand == "-"]
  genes$p3 <- genes$end
  genes$p3[genes$strand == "-"] <- genes$start[genes$strand == "-"]
  genes$size <- genes$end - genes$start
  cat("done (", nrow(genes), "genes)\n")

  cat("parsing exon information ... ")
  exons <- do.call(rbind, lapply(1:nrow(x), function(i) {
    df <- do.call(cbind, strsplit(as.character(x[i, c(10, 11)]), ","))
    cbind(df, rep(as.character(x[i, 1]), nrow(df)))
  }))
  exons <- data.frame(gene = as.character(exons[, 3]), start = as.numeric(exons[, 1]), end = as.numeric(exons[, 2]), stringsAsFactors = F)
  exons$chr <- genes$chr[match(exons$gene, genes$name)]
  exons <- exons[!duplicated(paste(exons[, 1], exons[, 2], exons[, 3])), ]
  genes <- genes[order(genes$size, decreasing = T), ]
  genes <- genes[!duplicated(genes$name), ]
  cat("done\n")

  cat("reading in", length(bam.files), "bam files ... ")
  cdl <- parallel::mclapply(bam.files, velocyto.R:::t.annotate.bam.reads,
                            genes = genes, exons = exons, margin = 1, exon.margin = 1,
                            mc.cores = n.cores, mc.preschedule = F)
  cat("done\n")

  cat("estimating gene counts ... ")
  edl <- parallel::mclapply(cdl, velocyto.R:::t.get.estimates2, genes = genes, mc.cores = n.cores, mc.preschedule = F)
  cat("done\n")

  cat("adjusting gene annotation based on expressed regions ... ")
  tl <- Matrix::colSums(do.call(rbind, parallel::mclapply(cdl, function(x) {
    ect <- table(c(x$exonstart, x$exonend))
    fect <- rep(0, nrow(exons))
    fect[as.integer(names(ect))] <- ect
    fect
  }, mc.cores = n.cores, mc.preschedule = F)))
  expr.exons <- exons[tl > min.exon.count, ]
  expr.lstat <- velocyto.R:::lengthstats2(1e+08, genes = genes, exons = expr.exons)
  df <- data.frame(il = log10(expr.lstat[, 2] + 1), el = log10(expr.lstat[, 3] + 1))
  rownames(df) <- rownames(expr.lstat)
  df$nex <- as.integer(table(expr.exons$gene)[rownames(df)])
  df$nex[is.na(df$nex)] <- 0
  cat("done\n")

  emat <- do.call(cbind, lapply(edl, function(d) {
    (d[, "exon"])
  }))
  emat[!is.finite(emat)] <- 0
  emat <- methods::as(emat, "dgCMatrix")
  smat <- do.call(cbind, lapply(edl, function(d) {
    (d[, "span"])
  }))
  smat[!is.finite(smat)] <- 0
  smat <- methods::as(smat, "dgCMatrix")
  iomat <- do.call(cbind, lapply(edl, function(d) {
    (d[, "introno"])
  }))
  iomat[!is.finite(iomat)] <- 0
  iomat <- methods::as(iomat, "dgCMatrix")
  return(list(emat = emat, iomat = iomat, smat = smat, base.df = df,
              exons = exons, genes = genes, expr.lstat = expr.lstat))

}

####################################################### READ DATA
### get bam files, filter out smart-SEQ2 samples
# bam
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Mapped/STAR_mm10"
bam_files <-
  list.files(path = bam_path, pattern = "*.bam$", full.names = T, recursive = F) %>%
  .[!str_detect(., "liver|fibroblast")] %>%
  magrittr::set_names(., str_replace(., ".*\\/(.*).SE.genome.Aligned.sortedByCoord.out.bam", "\\1"))

# GEO table with samples which are sequenced by SMART-SEQ2
geo_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Documentation/Deng_2014_Science_GSE45719.runInfo.xml"
geo_df <-
  XML::xmlToDataFrame(geo_path, stringsAsFactors = F) %>%
  as.tibble(.) %>%
  dplyr::mutate(Ids = str_extract(Ids, "SRS\\d{6}")) %>%
  dplyr::select(SRA_Sample_s = Ids, Description) %>%
  dplyr::filter(str_detect(Description, "fibroblast|split|smartseq2")) %>%
  dplyr::arrange(Description)

# runInfo table with Run_s info
runinfo_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Documentation/Deng_2014_Science_GSE45719.runInfo.txt"
runinfo_df <-
  readr::read_delim(file = runinfo_path, delim = "\t") %>%
  dplyr::select(Run_s, SRA_Sample_s) %>%
  dplyr::filter(SRA_Sample_s %in% geo_df$SRA_Sample_s)

# rename script with Run_s - sample name relation
rename_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Raw/Links/rename_links.sh"
rename_df <-
  read_delim(file = rename_path, delim = " ", col_names = c("command", "run", "sample_name")) %>%
  dplyr::mutate(Run_s = str_replace(run, ".fastq.gz", ""),
                sample_name = str_replace(sample_name, ".SE.txt.gz", "")) %>%
  dplyr::filter(Run_s %in% runinfo_df$Run_s)

# sample table with sample names and colors, filter out smart-seq2, order by stage
sample_df <-
  tibble(sample_name = names(bam_files),
         sample_path = bam_files) %>%
  dplyr::filter(!(sample_name %in% rename_df$sample_name)) %>%
  dplyr::mutate(stage = str_replace_all(sample_name, "_r.*|s_", ""),
                stage = str_replace(stage, "_", " "),
                color = gg_color_hue(n = length(unique(stage)), times = c(table(stage)))) %>%
  dplyr::right_join(.,
                    tibble(stage = c("MII oocyte", "zygote", "early 2C_blastomere" , "mid 2C_blastomere", "late 2C_blastomere",
                                     "2C blastomere", "4C blastomere", "8C blastomere", "16C blastomere",
                                     "early blastocyst", "mid blastocyst", "late blastocyst")),
                    by = "stage") %>%
  dplyr::filter(stage != "2C blastomere") %T>%
  readr::write_delim(., path = "Deng_2014.oocyte_no2C.sampleTable.tsv", delim = "\t")

### read bam files
# filter bam files
bam_files <- bam_files[names(bam_files) %in% sample_df$sample_name]

# refFlat annotation
annot_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/refFlat.mm10.20171127.txt.gz"

# parse gene annotation, annotate bam file reads
dat <- read.smartseq2.bams.modified(bam.files = bam_files, annotation.file = annot_path, n.cores = 20)
saveRDS(dat, file = "Deng_2014.oocyte_no2C.readsmartseq2dat.rds")