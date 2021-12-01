library(plyr)
library(stringr)
library(data.table)
library(GenomicAlignments)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing")

GffToGRanges = function(gff, filter=NULL){
  gff <- gtf
  
  if(ncol(gff) != 9)
    stop("Number of columns does not match gff format")
  
  if(any(gff[, 5] < gff[, 4])){
    warning("gff file contains ranges with negative widths...")
    gff = gff[gff[, 5] > gff[, 4], ]
  }
  
  if(!is.null(filter)){
    if(filter %in% gff[,3]){
      gff = gff[gff[, 3] == filter, ]
    }else{
      stop("The given feature is not present in the gff file")
    }
  }
  
  s <- strsplit(as.character(gff$V9), split = ";")
  z <- sapply(s, length)
  a <- split(s, z)
  gff <- gff[order(z), ]
  
  l <- lapply(a, function(x){
    d = sub('^ ','', unlist(x, use.names=F))
    d = sub('^.+? ','', d)
    m = matrix(d, ncol = length(x[[1]]), byrow=T)
    colnames(m) = sub(' .+$', '', sub('^ ', '', x[[1]]))
    m})
  
  ids <- rbind.fill(lapply(l, data.frame))
  ids <- ids[, c("gene_id", "transcript_id", "exon_id")]
  gff$V7[!gff$V7 %in% c('+', '-')] <- '*'
  granges <- GRanges(seqnames = gff[, 1],
                     IRanges(gff[, 4], gff[, 5]),
                     strand = gff[, 7],
                     frame = gff[, 8],
                     feature.type = gff[, 3],
                     .id = 1:nrow(gff))
  
  values(granges) <- cbind(values(granges), DataFrame(ids)[granges$.id, ])
  values(granges)$.id <- NULL
  return(granges)
}

# reads in the ensembl annotation, a bit of fixing some values
gtf <- read.table("/common/WORK/fhorvat/reference/mm10/Ensembl_GRCm38.86.20161128.gtf.gz", header = FALSE, sep = '\t', stringsAsFactors = F)
gtf <- gtf[!str_detect(gtf[ ,1], "GL|JH"), ]
gtf[gtf[, 1] == "MT", 1] <- "M"
gtf[, 1] <- paste0("chr", gtf[, 1])

# make GRanges from gtf (filter only exons)
gtf <- GffToGRanges(gtf, 'exon')

# get unique gene_id-transcript_id combinations
gids <- unique(values(gtf)[c('gene_id', 'transcript_id')])

# splits exon ranges on transcripts
gtf.trans <- split(gtf, gtf$transcript_id)

# orders transcripts based on number of exons in transcripts
gtf.trans <- gtf.trans[order(-elementNROWS(gtf.trans))] 

# keeps only first transcript of each gene (the one with most exons)
gtf.trans <- gtf.trans[!duplicated(gids$gene_id[match(names(gtf.trans), gids$transcript_id)])]

# gets only transcripts with 1 exon
gtf.trans.single <- gtf.trans[elementNROWS(gtf.trans) == 1]
gtf.trans.single <- unlist(gtf.trans.single)
gtf.trans.single$ex.num <- 1
gtf.trans.single$ex.tot <- 1
gtf.trans.single$gene_id <- NULL
gtf.trans.single$exon_id <- NULL
gtf.trans.single$feature.type <- NULL
gtf.trans.single$frame <- NULL
gtf.trans.single <- split(gtf.trans.single, gtf.trans.single$transcript_id)

# removes transcripts with more than 1 exon, reduces exons (merges them if they are closer than 100 bp)
gtf.trans <- gtf.trans[elementNROWS(gtf.trans) > 1]
gtf.trans <- reduce(gtf.trans, min.gapwidth = 100)
gtf.trans <- unlist(gtf.trans)
gtf.trans$transcript_id <- names(gtf.trans)

# counts exons in each transcript and enumerates them based on strand
d.val <- data.table(as.data.frame(values(gtf.trans)))
d.val$strand <- as.character(strand(gtf.trans))
d.val[d.val$strand == '+' , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == '+']]
d.val[d.val$strand == '-' , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == '-']]
gtf.trans$ex.num <- d.val$IDX
gtf.trans$ex.tot <- d.val$COUNT
gtf.trans <- split(gtf.trans, as.character(gtf.trans$transcript_id))

# gets exons
gtf.ex <- unlist(gtf.trans)

# constructs the splicing coordinates
splice <- unlist(gtf.trans)
names(splice) <- NULL

# splice donor
splice.don <- GenomicRanges::resize(splice, fix = "end", width = 1)
# splice.don <- c(GenomicRanges::resize(splice.don[strand(splice.don) == "+"], width = 3, fix = "start"), 
#                 GenomicRanges::resize(splice.don[strand(splice.don) == "-"], width = 3, fix = "end"))
# splice.don <- GenomicRanges::shift(splice.don, -1)
splice.don <- splice.don[splice.don$ex.num != splice.don$ex.tot] #removes the last exon in the transcript from splice donors
splice.don <- unique(splice.don)
splice.don <- splice.don[countOverlaps(splice.don, splice.don) == 1] #removes splice donors which overlap each other
write.table(splice.don, file = "Ensembl_GRCm38.86.20161128_splice_donor.txt", sep = "\t", row.names = F)

# splice acceptor
splice.acc <- GenomicRanges::resize(splice, fix = "start", width = 1)
# splice.acc <- c(GenomicRanges::resize(splice.acc[strand(splice.acc) == "+"], width = 3, fix = "end"), 
#                 GenomicRanges::resize(splice.acc[strand(splice.acc) == "-"], width = 3, fix = "start"))
# splice.acc <- GenomicRanges::shift(splice.acc, 1)
splice.acc <- splice.acc[splice.acc$ex.num != 1] #removes the first exon in the transcript from splice acceptors
splice.acc <- unique(splice.acc)
splice.acc <- splice.acc[countOverlaps(splice.acc, splice.acc) == 1]
names(splice.acc) <- NULL
write.table(splice.acc, file = "Ensembl_GRCm38.86.20161128_splice_acceptor.txt", sep = "\t", row.names = F)

# joins transcripts with one and with more exons, gets ranges, finds introns
gtf.genes <- c(gtf.trans, gtf.trans.single)
write.table(unname(unlist(gtf.genes)), file = "Ensembl_GRCm38.86.20161128_exons.txt", sep = "\t", row.names = F)
gtf.range <- unlist(range(gtf.genes))
gtf.int <- GenomicRanges::setdiff(gtf.range, gtf.ex)

