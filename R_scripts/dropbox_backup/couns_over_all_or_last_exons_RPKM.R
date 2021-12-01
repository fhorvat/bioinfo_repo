library("Rsamtools")
library("GenomicFeatures")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("GenomicRanges")
library("BiocParallel")
library("GenomicAlignments")
library("DESeq2")

countingAllorLastExon <- function(ebg_features, bam_filenames, out_name){
  bamfiles <- BamFileList(bam_filenames, yieldSize = 2000000)
  sample_table <- DataFrame(sample_names = gsub("_Aligned.sortedByCoord.out.bam|mapping/STAR/", "", bam_filenames), 
                            row.names = gsub("_Aligned.sortedByCoord.out.bam|mapping/STAR/", "", bam_filenames))
  
  # counting 
  register(MulticoreParam())
  exons_se <- summarizeOverlaps(features = ebg_features, 
                                reads = bamfiles, 
                                mode = "Union", 
                                singleEnd = TRUE, 
                                ignore.strand = TRUE)
  colData(exons_se) <- sample_table
  if (out_name == "last"){
    rownames(assay(exons_se)) <- names(all_exons_ebg)
  }
  
  # counts table
  exons_counts <- data.frame(assay(exons_se))
  write.csv(exons_counts, paste0(out_name, "_exons_counts.csv"))
  
  # rpkm / fpkm table
  exons_dds <- DESeqDataSet(exons_se, ~ 1)
  exons_rpkm <- fpkm(exons_dds,  robust = FALSE)
  write.csv(exons_rpkm, paste0(out_name, "_exons_rpkm.csv"))
  
  return("done")
}

countingAllorLastExonByPath <- function(path, exon_all_or_last){
  setwd(path) 
  filenames <- list.files(path = "mapping/STAR", pattern = "\\.bam$", full.names = T)
  if (exon_all_or_last == "all"){
    countingAllorLastExon(ebg_features = all_exons_ebg, bam_filenames = filenames, out_name = "all")
  } else{
    if (exon_all_or_last == "last"){
      countingAllorLastExon(ebg_features = last_exons_ebg, bam_filenames = filenames, out_name = "last")  
    } else{
      stop("must be all or last")
    }
  }
  return("done")
}

# making ebg objects for counting (all exon or only last exon in gene)
all_exons_ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm9.knownGene, by = "gene")
last_exons_ebg <- lapply(X = 1:length(all_exons_ebg), function(X) all_exons_ebg[[X]][length(all_exons_ebg[[X]])])
last_exons_ebg <- GRangesList(last_exons_ebg)
names(last_exons_ebg) <- names(all_exons_ebg)

# function call
countingAllorLastExonByPath(path = "/common/WORK/fhorvat/Projekti/Svoboda/Other/MousePolyA/MII/Hamazaki_2015",
                            exon_all_or_last = "last")
