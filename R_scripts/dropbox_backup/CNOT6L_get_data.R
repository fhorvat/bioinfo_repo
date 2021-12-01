### INFO: Takes pair end reads and getws intron/exon splice ratio
### DATE: 08. 03. 2017
### AUTHOR: Vedran Franke, Filip Horvat
rm(list=ls()); gc()

# WORKING DIR 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads")

# LIBRARIES
lib.path <- '/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/vfranke_documentation'
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'ScanLib.R'))
library(doMC)
library(GenomicAlignments)
library(genomation)
library(data.table)
library(plyr)
library(readr)
library(stringr)
library(tibble)
library(magrittr)

# PATHS
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads"
r.outpath <- file.path(outpath, 'R_objects');
# dir.create(r.outpath, showWarnings=FALSE)

# SCRIPT PARAMETERS
registerDoMC(21)

################################################################## experiment table
# CNOT6L experiment table
sample_table <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_sample_list_11919R_2015_10_29.csv", col_names = T) %>%
  dplyr::select(ID, stage = `Time Course`, treatment = `Treatment/Control`) %>%
  dplyr::mutate(name = str_c(ID, stage, treatment, sep = "_")) 

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
  dplyr::mutate(ID = str_replace_all(sample_path, "^/.*/|_.*", "")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(library_size = as.integer(read_delim(logs_path, delim = "\t", col_names = F)[8, 2])) %>% 
  dplyr::select(ID, library_size, sample_path)

# join together, select only 1C knock-out (for now)
sample_table <- 
  dplyr::left_join(sample_table, library_size_df, by = "ID") %>% 
  dplyr::filter(treatment == "WT")

sample_path <- sample_table$sample_path
################################################################## 
# repeatMasker 
reps <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# reads in the ensembl annotation
gtf = read.table("/common/WORK/fhorvat/reference/mouse/mm10/Ensembl_GRCm38.86.20161128.gtf.gz", header=FALSE, sep='\t', stringsAsFactors = F)
gtf[, 1] = paste0("chr", gtf[, 1])
gtf = gtf[!str_detect(gtf[,1],'NT'),]
gtf[gtf[,1] == 'chrMT',1] = 'chrM'
gtf = GffToGRanges(gtf, 'exon')

gids = unique(values(gtf)[c('gene_id','transcript_id')])
gtf.trans = split(gtf, gtf$transcript_id)
gtf.trans = gtf.trans[order(-elementNROWS(gtf.trans))]
gtf.trans = gtf.trans[!duplicated(gids$gene_id[match(names(gtf.trans), gids$transcript_id)])]

gtf.trans.single = gtf.trans[elementNROWS(gtf.trans) == 1]
gtf.trans.single = unlist(gtf.trans.single)
values(gtf.trans.single) = gtf.trans.single$transcript_id
names(mcols(gtf.trans.single)) = "transcript_id"
gtf.trans.single$ex.num = 1
gtf.trans.single$ex.tot = 1
gtf.trans.single = split(gtf.trans.single, gtf.trans.single$transcript_id)

gtf.trans = gtf.trans[elementNROWS(gtf.trans) > 1]
gtf.trans = reduce(gtf.trans, min.gapwidth=100)
gtf.trans = unlist(gtf.trans)
gtf.trans$transcript_id = names(gtf.trans)
d.val = data.table(as.data.frame(values(gtf.trans)))
d.val$strand = as.character(strand(gtf.trans))
d.val[d.val$strand == '+' , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == '+']]
d.val[d.val$strand == '-' , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == '-']]
gtf.trans$ex.num = d.val$IDX
gtf.trans$ex.tot = d.val$COUNT
gtf.trans = split(gtf.trans, as.character(gtf.trans$transcript_id))
gtf.ex = unlist(gtf.trans)

# constructs the splicing coordinates
splice = unlist(gtf.trans)
splice.don <- GenomicRanges::resize(splice, fix = "end", width = 1)
splice.don = splice.don[splice.don$ex.num != splice.don$ex.tot]
splice.don = resize(splice.don, width=10, fix='start')
splice.don = unique(splice.don)
splice.don = splice.don[countOverlaps(splice.don, splice.don) == 1]

splice.acc <- GenomicRanges::resize(splice, fix = "start", width = 1)
splice.acc = splice.acc[splice.acc$ex.num != 1]	
splice.acc = resize(splice.acc, width=10, fix='end')
splice.acc = unique(splice.acc)
splice.acc = splice.acc[countOverlaps(splice.acc, splice.acc) == 1]

gtf.genes = c(gtf.trans, gtf.trans.single)
gtf.range = unlist(range(gtf.genes))
gtf.int = setdiff(gtf.range, gtf.ex)

# loops through the file and calculates the statistics
for(file in sample_path){
  
  name = BamName(file)
  print(name)
  chrs = chrFinder(file)
  chrs = chrs[chrs$chr != 'chrM',]
  
  cnts = foreach(chr = chrs$chr)%dopar%{
    
    print(chr)
    
    # reads in the reads
    cat('Reading files...\n')
    which.ranges = GRanges(chr, IRanges(1, chrs$chrlen[chrs$chr == chr]))
    if(str_detect(name, 'SE')){
      bam = readGAlignments(file, param=ScanBamParam(which=which.ranges, tag='NH'), use.names=TRUE)
      bams = bam
    }else{
      bam = suppressWarnings(readGAlignmentPairs(file, param=ScanBamParam(which=which.ranges), use.names=TRUE))
      bams = readGAlignments(file, param=ScanBamParam(which=which.ranges, tag='NH'), use.names=TRUE)
    }
    
    gbam = grglist(bam)
    tab = unique(data.table(names=names(bams), V1=values(bams)$NH))
    tab = tab[tab$names %in% names(bam)]
    tab = tab[tab$V1==1]
    rm(bam, bams);gc()
    
    cat('Selecting reads...\n')
    co.ex = countOverlaps(gbam, gtf.ex, ignore.strand=TRUE)
    co.in = countOverlaps(gbam, gtf.int, ignore.strand=TRUE, minoverlap=15)
    co.ind = which(co.ex > 1 | (co.ex >= 1 & co.in >= 1))
    
    forep = data.table(as.matrix(findOverlaps(gbam, reps, ignore.strand=TRUE)))
    
    fot = data.table(as.matrix(findOverlaps(gbam, gtf.range, ignore.strand=TRUE)))
    fot$id = names(gbam)[fot$queryHits]
    fon = fot[,length(unique(subjectHits)), by=id]
    fon = fon[fon$V1 == 1]
    fot = fot[fot$id %in% fon$id]

    fot$ens.trans.id = names(gtf.range)[fot$subjectHits]
    fot$dist = max(start(gbam)[fot$queryHits] - min(end(gbam)[fot$queryHits]))
    fot$dist.2 = max(end(gbam)[fot$queryHits]) - min(start(gbam)[fot$queryHits])
    fot$reps = fot$queryHits %in% forep$queryHits
    fot = fot[fot$queryHits %in% co.ind]
    fot = fot[fot$id %in% tab$names]
    fot = fot[!fot$reps]
    fot$ens.trans.id = NULL
    fot$subjectHits = NULL
    fot = unique(fot)
    fot$ex = fot$queryHits %in% which(co.ex>1)
    fot$int = fot$queryHits %in% which(co.in>0)
    fot$ex[fot$int > 0] = FALSE

    cat('Counting ex - int...\n')
    fo.don = data.table(as.matrix(findOverlaps(splice.don, gbam, type='within')))
    fo.don$names = names(gbam)[fo.don$subjectHits]
    fo.don$ens.trans.id = splice.don$transcript_id[fo.don$queryHits]
    fo.acc = data.table(as.matrix(findOverlaps(splice.acc, gbam, type='within')))
    fo.acc$names = names(gbam)[fo.acc$subjectHits]
    fo.acc$ens.trans.id = splice.acc$transcript_id[fo.acc$queryHits]
    
    cat('Marking reads ex - int...\n')
    tab$don = tab$names %in% fo.don$names
    tab$acc = tab$names %in% fo.acc$names
    tab = tab[tab$don | tab$acc]
    
    cat('Counting reads ex - int...\n')
    fo.don.co = fo.don[,length(subjectHits), by=ens.trans.id]
    fo.acc.co = fo.acc[,length(subjectHits), by=ens.trans.id]
    m = merge(fo.don.co, fo.acc.co, by='ens.trans.id', all=T)
    m[is.na(m)] = 0
    setnames(m,2:3, c('don','acc'))
    
    return(list(tab=tab,m = m,dist=fot))
  }
  
  tab = lapply(cnts, '[[', 'tab')
  d.reads = rbindlist(tab)
  d.reads = d.reads[order(-rowSums(d.reads[,c('don','acc'), with=F]))]
  d.reads = d.reads[!duplicated(d.reads$names)]
  
  dist = rbindlist(lapply(cnts, '[[', 'dist'))
  
  cat('Saving data...\n')
  l = list(d.reads, dist)
  save(l, file = file.path(r.outpath, paste(name,'SpliceDonAcc.RData',sep='.')))
  
}

