library(dplyr)
library(readr)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tibble)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(BiocParallel)
library(doMC)
library(genomation)

### replicating figure 6D from EMBO-J paper - comparison of unspliced/spliced read pair ratios per cell stage
# UNSPLICED 
# - read counts where one end maps to intron/exon junction or entirely in intron and the
#   other end maps to the adjacent exon were labeled as ‘unspliced’
# SPLICED 
# - one end maps either to the splice site and covers two adjacent exons or 
#   with each end mapping to separate, adjacent exons

# some code from vfranke:
# /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/SplicingAnalysis/SplicingAnalysis_Scripts

################################################################## functions
# modified vfranke function, takes gff and outputs GRanges (with optional filtering)
GffToGRanges <- function(gff, filter = NULL){
  
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(tibble)
  
  # gff has to have 9 columns
  if(ncol(gff) != 9)
    stop("Number of columns does not match gff format")
  
  # width of elements have to be positive
  if(any(as.integer(gff$X5) < as.integer(gff$X4))){
    warning("gff file contains ranges with negative widths...")
    gff <- gff[as.integer(gff$X5) > as.integer(gff$X4), ]
  }
  
  # filter - for example: exon/transcript/gene/CDS...
  if(!is.null(filter)){
    if(filter %in% gff$X3){
      gff = gff[gff$X3 == filter, ]
    }else{
      stop("The given feature is not present in the gff file")
    }
  }
  
  ### get exon/transcript/gene ID for the 9th column of gff table
  # split stings in 9th column, split list again by length (not all genes have all info in that column)
  s <- str_split(string = gff$X9, pattern = ";")
  z <- sapply(s, length)
  a <- split(s, z)
  gff <- gff[order(z), ]
  
  # create tables from 9th column with different number of columns in a list
  l <- lapply(a, function(x){
    d <- str_trim(str_replace(unlist(x, use.names = F), "^.+? ", "")) %>% 
      str_replace_all(., "\"", "")
    m <- matrix(d, ncol = length(x[[1]]), byrow = T) %>% 
      set_colnames(str_replace(str_trim(x[[1]]), " .+$", "")) %>% 
      as_tibble()
    return(m)
  })
  
  # bind all tables from list in one with filling empty data with NA 
  ids <- 
    dplyr::bind_rows(l) %>% 
    dplyr::select(gene_id, transcript_id, exon_id)
  
  # set all strands which are not +/- to *
  gff$X7[!gff$X7 %in% c('+', '-')] <- "*"
  
  # create and output GRanges
  granges <- GRanges(seqnames = gff$X1,
                     IRanges(as.integer(gff$X4), as.integer(gff$X5)),
                     strand = gff$X7,
                     frame = gff$X8,
                     feature.type = gff$X3,
                     .id = 1:nrow(gff))
  
  values(granges) <- cbind(values(granges), DataFrame(ids)[granges$.id, ])
  values(granges)$.id <- NULL
  
  return(granges)
}

BamName <-  function(bam.file){
  sub('.bam', '', basename(bam.file))
}

chrFinder <- function(bam.path, filter = FALSE, output = "data.frame"){
  
  s <- scanBamHeader(bam.path)
  st <- s[[1]]$text
  st <- do.call(rbind,st[names(st) == "@SQ"])
  st[, 1] <- str_replace(st[, 1], "SN:", "")
  st[, 2] <- str_replace(st[, 2], "LN:", "")
  
  if(filter == TRUE){
    st <- st[!str_detect(st[, 1], "random")]
  }
  
  if(output == "data.frame"){
    vst <- data.frame(chr = st[, 1], chrlen = as.numeric(st[, 2]))
  }else{
    vst <- as.numeric(st[, 2])
    names(vst) <- st[, 1]
  }
  return(vst)
}

################################################################## parameters 
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads")

#  register the multicore parallel backend with the foreach package 
registerDoMC(21)

# set paths
r.outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/R_objects"

################################################################## getting genes, exons, introns, promoters and intergenic regions from knownGene tables
# read ENSEMBL gene table
ensembl_genes <- read_delim(file = "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl_GRCm38.86.20161128.gtf.gz", 
                            delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# filter and clean, make GRanges from gtf (filter only exons)
gtf <- 
  ensembl_genes %>% 
  dplyr::filter(!str_detect(X1, pattern = "GL|JH")) %>% 
  mutate(X1 = replace(X1, X1 == "MT", "M"), 
         X1 = paste0("chr", X1)) %>% 
  GffToGRanges(., filter = "exon")

################################################################## experiment table, input reads
# CNOT6L experiment table
sample_table <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_sample_list_11919R_2015_10_29.csv", col_names = T) %>%
  dplyr::select(ID, stage = `Time Course`, treatment = `Treatment/Control`) %>%
  mutate(name = str_c(ID, stage, treatment, sep = "_")) %>% 
  dplyr::filter(stage == "1C", treatment == "KO")

# CNOT6L library size
library_size_df <- 
  tibble(sample_path = list.files("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                  pattern = "*.bam$",
                                  recursive = T, 
                                  full.names = T), 
         logs_path = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                pattern = "*Log.final.out", 
                                recursive = T, 
                                full.names = T)) %>%
  mutate(ID = str_replace_all(sample_path, "^/.*/|_.*", "")) %>%
  rowwise() %>%
  mutate(library_size = as.integer(read_delim(logs_path, delim = "\t", col_names = F)[8, 2])) %>% 
  dplyr::select(ID, library_size, sample_path)

# join together, select only 1C knock-out (for now)
sample_table <- 
  left_join(sample_table, library_size_df, by = "ID") %>% 
  dplyr::filter(stage == "1C", treatment == "KO")

################################################################## 
# get unique gene_id-transcript_id combinations
gids <- unique(values(gtf)[c('gene_id', 'transcript_id')])

# splits exon ranges on transcripts
gtf.trans <- split(gtf, gtf$transcript_id)

# orders transcripts based on number of exons in transcripts
gtf.trans <- gtf.trans[order(-elementNROWS(gtf.trans))] 

# keeps only first transcript of each gene (the one with most exons)
gtf.trans <- gtf.trans[!duplicated(gids$gene_id[match(names(gtf.trans), gids$transcript_id)])]

# gets only transcripts with 1 exon
gtf.trans.single <- 
  gtf.trans[elementNROWS(gtf.trans) == 1] %>% 
  unlist(.)
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

################################################################## splice cooridnates 
# constructs the splicing coordinates
splice <- unlist(gtf.trans)

# splice donor
splice.don <- GenomicRanges::resize(splice, fix = "end", width = 1)
splice.don <- splice.don[splice.don$ex.num != splice.don$ex.tot] # removes the last exon in the transcript from splice donors
splice.don <- resize(splice.don, width = 10, fix = "start")
splice.don <- unique(splice.don)
splice.don <- splice.don[countOverlaps(splice.don, splice.don) == 1] # removes splice donors which overlap each other

# splice acceptor
splice.acc <- GenomicRanges::resize(splice, fix = "start", width = 1)
splice.acc <- splice.acc[splice.acc$ex.num != 1] # removes the first exon in the transcript from splice acceptors
splice.acc <- resize(splice.acc, width = 10, fix = "end")
splice.acc <- unique(splice.acc)
splice.acc <- splice.acc[countOverlaps(splice.acc, splice.acc) == 1]

# joins transcripts with one and with more exons, gets ranges, finds introns
gtf.genes <- c(gtf.trans, gtf.trans.single)
gtf.range <- unlist(range(gtf.genes))
gtf.int <- GenomicRanges::setdiff(gtf.range, gtf.ex)

################################################################## 

# goes through files in a loop
for(n in nrow(sample_table)){
  
  # loops through the file and calculates the statistics
  file <- as.character(sample_table[n, "sample_path"])
  name <- BamName(file)
  print(name)
  chrs <- chrFinder(file) %>% 
    dplyr::filter(chr != "chrM")
  
  # does parallel computing 
  cnts <- foreach(chr = chrs$chr) %dopar% {
    
    print(chr)
    
    # reads in the reads
    cat('Reading files...\n')
    which.ranges <- GRanges(chr, IRanges(1, chrs$chrlen[chrs$chr == chr]))
    
    if(str_detect(name, "SE")){
      bam <- readGAlignments(file, param = ScanBamParam(which = which.ranges, tag = "NH"), use.names = TRUE)
      bams <- bam
    }else{
      bam <- suppressWarnings(readGAlignmentPairs(file, param = ScanBamParam(which = which.ranges), use.names = TRUE))
      bams <- readGAlignments(file, param = ScanBamParam(which = which.ranges, tag = "NH"), use.names = TRUE)
    }
    
    # split reads to list (if is paired-end than each fragment is split to two ranges)
    gbam <- grglist(bam)
    
    # get names of uniquely mapped reads (NH = 1)
    tab <- unique(data.table(names = names(bams), V1 = values(bams)$NH))
    tab <- tab[tab$names %in% names(bam)]
    tab <- tab[tab$V1 == 1]
    rm(bam, bams); gc()
    
    # counts overlaps between reads and exons/introns
    # gets indices of reads which overlap either:
    # - one exon
    # - one exon and one intron
    cat('Selecting reads...\n')
    co.ex <- countOverlaps(gbam, gtf.ex, ignore.strand = TRUE)
    co.in <- countOverlaps(gbam, gtf.int, ignore.strand = TRUE, minoverlap = 15)
    co.ind <- which(co.ex > 1 | (co.ex >= 1 & co.in >= 1))
    
    # # finds overlaps betweeen reads and repeatMasker
    # forep <- data.table(as.matrix(findOverlaps(gbam, reps, ignore.strand = TRUE)))
    
    # finds overlaps between reads and transcript
    fot <- data.table(as.matrix(findOverlaps(gbam, gtf.range, ignore.strand = TRUE)))
    fot$id <- names(gbam)[fot$queryHits]
    
    # counts how many transcripts reads overlap, take only those reads which overlap one range
    fon <- fot[, length(unique(subjectHits)), by = id]
    fon <- fon[fon$V1 == 1]
    fot <- fot[fot$id %in% fon$id]
    
    # add ID of transcript to reads which are mapped to them
    fot$ens.trans.id <- names(gtf.range)[fot$subjectHits]
    
    # adds T/F values if gene is on the chosen gene list
    # for(n in names(l.genes))
    #   fot[[n]] <- fot$ens.trans.id %in% l.genes[[n]]$ens.trans.id
    
    # calculates distance between two reads on one pair-end fragment  
    fot$dist <- max(start(gbam)[fot$queryHits] - min(end(gbam)[fot$queryHits]))
    
    # calculates length of the fragment
    fot$dist.2 <- max(end(gbam)[fot$queryHits]) - min(start(gbam)[fot$queryHits])
    
    # # adds T/F values if gene is overlaping repeatMasker
    # fot$reps <- fot$queryHits %in% forep$queryHits
    
    # filtering, keep reads: 
    # - which overlap one exon or one exon and one intron
    # - which are uniquely mapped (no multi-mappers)
    # - which don't overlap ranges from reapeatMasker
    # - which overlap only one transcript
    fot <- fot[fot$queryHits %in% co.ind]
    fot <- fot[fot$id %in% tab$names] 
    # fot <- fot[!fot$reps]
    # fot <- fot[fot$queryHits %in% fon$queryHits]
    
    # get unique combinations of reads-transcripts
    fot$ens.trans.id <- NULL
    fot$subjectHits <- NULL
    fot <- unique(fot)
    
    # add T/F for reads overlaping exons/introns, adds F exons if reads overlap introns
    fot$ex <- fot$queryHits %in% which(co.ex > 1)
    fot$int <- fot$queryHits %in% which(co.in > 0)
    fot$ex[fot$int > 0] <- FALSE
    # if(name == 's_MII.WE' & sum(fot$sel.samp.mii) != 0)
    # stop(paste(name, chr, 'smth went wrong'))
    
    cat('Counting ex - int...\n')
    # finding overlaps between splice donors and reads
    fo.don <- data.table(as.matrix(findOverlaps(splice.don, gbam, type = "within")))
    fo.don$names <- names(gbam)[fo.don$subjectHits]
    fo.don$ens.trans.id <- splice.don$transcript_id[fo.don$queryHits]
    
    # finding overlaps between splice acceptors and reads
    fo.acc <- data.table(as.matrix(findOverlaps(splice.acc, gbam, type = "within")))
    fo.acc$names <- names(gbam)[fo.acc$subjectHits]
    fo.acc$ens.trans.id <- splice.acc$transcript_id[fo.acc$queryHits]
    
    cat('Marking reads ex - int...\n')
    # marking if reads are mapped to splice donor/acceptor, remove ones which are not mapped to any
    tab$don <- tab$names %in% fo.don$names
    tab$acc <- tab$names %in% fo.acc$names
    tab <- tab[tab$don | tab$acc]
    
    cat('Counting reads ex - int...\n')
    # counting how many reads are mapped to splice donor/acceptor of particular transcripts
    fo.don.co <- fo.don[, length(subjectHits), by = ens.trans.id]
    fo.acc.co <- fo.acc[, length(subjectHits), by = ens.trans.id]
    m <- merge(fo.don.co, fo.acc.co, by = 'ens.trans.id', all = T)
    m[is.na(m)] <- 0
    setnames(m, 2:3, c('don', 'acc'))
    
    return(list(tab = tab, m = m, dist = fot))
  }
  
  # extract reads which are mapped to donor or acceptor
  tab <- lapply(cnts, "[[", "tab")
  d.reads <- rbindlist(tab)
  d.reads <- d.reads[order(-rowSums(d.reads[, c('don','acc'), with = F]))]
  d.reads <- d.reads[!duplicated(d.reads$names)]
  
  # m <- rbindlist(lapply(cnts, '[[', 'm'))
  # setnames(m, 1:3, c('ens.trans.id','don','acc'))
  
  # extract data about reads (distances between two reads on fragment, fragment length, mapping to exon/intron)
  dist <- rbindlist(lapply(cnts, "[[", "dist"))
  
  # save data as R object
  cat('Saving data...\n')
  l <- list(d.reads, dist)
  save(l, file = file.path(r.outpath, paste(name, "SpliceDonAcc.RData", sep = ".")))
  
}

