#!/common/WORK/vfranke/bin/R/R-3.1.0/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/no_genes_sel")

# {1} LIBRARIES
lib.path='/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/vfranke_documentation'
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'ScanLib.R'))
library(data.table)
library(stringr)
library(doMC)
library(GenomicAlignments)
library(genomation)
library(readr)

#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS

#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES 

# {{{1}}} PATH VARIABLES
# file = '/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_MII.WE/s_MII.WE.bam'
# file = '/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_1cell.WE/s_1cell.WE.bam'
# file = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam'

sample_path <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_MII.WE/s_MII.WE.bam",
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE_DNAm/s_1cell.WE_DNAm.bam",
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam", 
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_Molura.WE/s_Molura.WE.bam", 
                           "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_Blast.WE/s_Blast.WE.bam"))

# outpath='/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/GenomeActivation/SplicingAnalysis'
outpath = '/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/no_genes_sel'

r.outpath = file.path(outpath, 'R_objects');
dir.create(r.outpath, showWarnings=FALSE)

mm9.gtf.path = '/common/WORK/fhorvat/reference/mouse/mm10/Ensembl_GRCm38.86.20161128.gtf.gz'

# ens.genes.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/25012012.EnsemblGenesAll_BioMart.form.txt'

genes.sel.path = '/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/vfranke_documentation/GeneCount_RepsRemoved_RPKM.txt'

# rep.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.RepeatMasker.RData'

#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
registerDoMC(21)

#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES


# {{3}} MAIN

# Assigner(rep.path, 'reps')

# repeatMasker table
reps <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# reads in the ensembl annotation
gtf = read.table(mm9.gtf.path, header=FALSE, sep='\t', stringsAsFactors = F)
gtf[, 1] = paste0("chr", gtf[, 1])
gtf = gtf[!str_detect(gtf[,1],'NT'),]
gtf[gtf[,1] == 'chrMT',1] = 'chrM'
gtf = GffToGRanges(gtf, 'exon')

gids = unique(values(gtf)[c('gene_id','transcript_id')])
gtf.trans = split(gtf, gtf$transcript_id)
# gtf.trans = gtf.trans[order(-sum(width(gtf.trans)))]
gtf.trans = gtf.trans[order(-elementNROWS(gtf.trans))]
gtf.trans = gtf.trans[!duplicated(gids$gene_id[match(names(gtf.trans), gids$transcript_id)])]

gtf.trans.single = gtf.trans[elementNROWS(gtf.trans) == 1]
gtf.trans.single = unlist(gtf.trans.single)
values(gtf.trans.single) = gtf.trans.single$transcript_id
names(mcols(gtf.trans.single)) = "transcript_id"
gtf.trans.single$ex.num = 1
gtf.trans.single$ex.tot = 1
gtf.trans.single = split(gtf.trans.single, gtf.trans.single$transcript_id)

# a = 'ENSMUST00000067468'
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

# # reads in the selected genes
# genes.sel = read.table(genes.sel.path, header=TRUE, dec=',', sep='\t')
# genes.sel = genes.sel[genes.sel$biotype == 'protein_coding',]
# 
# genes.sel.int = with(genes.sel, genes.sel[order(-genes.sel$s_1cell.WE.int.reps.rem),][1:1000,])
# genes.sel.rat = with(genes.sel, genes.sel[order(-log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1))),][1:1000,])
# 
# genes.sel.samp.rat = genes.sel[genes.sel$s_1cell.WE.int.reps.rem != 0,]
# genes.sel.samp.rat = with(genes.sel.samp.rat, genes.sel.samp.rat[order(-log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) /
#                                                                          log2((s_MII.WE.int.reps.rem+1)/(s_MII.WE.ex+1))),][1:1000,])
# genes.sel.mii.ex = subset(genes.sel, s_MII.WE.ex == 0 & s_1cell.WE.int.reps.rem != 0)
# genes.sel.mii.ex = genes.sel.mii.ex[order(-genes.sel.mii.ex$s_1cell.WE.int.reps.rem),]
# 
# genes.sel.mii.int = subset(genes.sel, s_MII.WE.int.reps.rem == 0 & s_1cell.WE.int.reps.rem != 0)
# genes.sel.mii.int = genes.sel.mii.int[order(-genes.sel.mii.int$s_1cell.WE.int.reps.rem),]
# 
# genes.sel.mii = subset(genes.sel, s_MII.WE.int.reps.rem == 0 & s_1cell.WE.int.reps.rem != 0 & s_MII.WE.ex == 0)
# 
# l.genes = list(genes.sel.int = genes.sel.int,
#                genes.sel.rat = genes.sel.rat,
#                genes.sel.samp.rat = genes.sel.samp.rat,
#                genes.sel.mii.ex = genes.sel.mii.ex,
#                genes.sel.mii.int = genes.sel.mii.int,
#                genes.sel.mii = genes.sel.mii)

# constructs the splicing coordinates
splice = unlist(gtf.trans)
# splice.don = c(resize(splice[strand(splice) == '+'], fix='end', width=1),
# 			   resize(splice[strand(splice) == '-'], fix='start', width=1))
splice.don <- GenomicRanges::resize(splice, fix = "end", width = 1)
splice.don = splice.don[splice.don$ex.num != splice.don$ex.tot]
splice.don = resize(splice.don, width=10, fix='start')
splice.don = unique(splice.don)
splice.don = splice.don[countOverlaps(splice.don, splice.don) == 1]

# splice.acc = c(resize(splice[strand(splice) == '+'], fix='start', width=1),
# 			   resize(splice[strand(splice) == '-'], fix='end', width=1))
splice.acc <- GenomicRanges::resize(splice, fix = "start", width = 1)
splice.acc = splice.acc[splice.acc$ex.num != 1]	
splice.acc = resize(splice.acc, width=10, fix='end')
splice.acc = unique(splice.acc)
splice.acc = splice.acc[countOverlaps(splice.acc, splice.acc) == 1]

gtf.genes = c(gtf.trans, gtf.trans.single)
gtf.range = unlist(range(gtf.genes))
gtf.int = setdiff(gtf.range, gtf.ex)

for(file in sample_path){
  
  # loops through the file and calculates the statistics
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
    # for(n in names(l.genes))
    # 	fot[[n]] = fot$ens.trans.id %in% l.genes[[n]]$ens.trans.id

    fot$dist = max(start(gbam)[fot$queryHits] - min(end(gbam)[fot$queryHits]))
    fot$dist.2 = max(end(gbam)[fot$queryHits]) - min(start(gbam)[fot$queryHits])
    fot$reps = fot$queryHits %in% forep$queryHits
    fot = fot[fot$queryHits %in% co.ind]
    fot = fot[fot$id %in% tab$names]
    fot = fot[!fot$reps]
    # fot = fot[fot$queryHits %in% fon$queryHits]
    fot$ens.trans.id = NULL
    fot$subjectHits = NULL
    fot = unique(fot)
    fot$ex = fot$queryHits %in% which(co.ex>1)
    fot$int = fot$queryHits %in% which(co.in>0)
    fot$ex[fot$int > 0] = FALSE
    # if(name == 's_MII.WE' & sum(fot$sel.samp.mii) != 0)
    # stop(paste(name, chr, 'smth went wrong'))
    
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
  
  # m = rbindlist(lapply(cnts, '[[', 'm'))
  # setnames(m, 1:3, c('ens.trans.id','don','acc'))
  
  dist = rbindlist(lapply(cnts, '[[', 'dist'))
  
  cat('Saving data...\n')
  l = list(d.reads, dist)
  save(l, file = file.path(r.outpath, paste(name,'SpliceDonAcc.RData',sep='.')))
  
}

cat('Thank you, and goodbye!...\n')
#/{2} CODE


# pdf(file.path(outpath, 'ExInt.Rat.pdf'))
# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_MII.WE.int.reps.rem+1)/(s_MII.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs MII WE'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_1cell.WE_DNAm.int.reps.rem+1)/(s_1cell.WE_DNAm.ex+1)), pch=20, cex=.75, main='1cell WE vs 1cell DRB'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_2cell.WE.int.reps.rem+1)/(s_2cell.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs 2cell WE'))

# with(genes.sel, plot(log2((s_1cell.PA.int.reps.rem + 1)/(s_1cell.PA.ex + 1)) ~ log2((s_MII.PA.int.reps.rem+1)/(s_MII.PA.ex+1)), pch=20, cex=.75, main='1cell PA vs 2cell PA'))
# dev.off()


# pdf(file.path(outpath, 'ExInt.Rat.1cell.rem.pdf'))
# gene.sel.MII = subset(genes.sel, s_MII.WE.int.reps.rem == 0 & s_MII.WE.ex ==0)
# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_MII.WE.int.reps.rem+1)/(s_MII.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs MII WE', col=rgb(0,0,0,alpha=.1)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_MII.WE.int.reps.rem+1)/(s_MII.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_1cell.WE_DNAm.int.reps.rem+1)/(s_1cell.WE_DNAm.ex+1)), pch=20, cex=.75, main='1cell WE vs 1cell DRB', col=rgb(0,0,0,alpha=.5)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_1cell.WE_DNAm.int.reps.rem+1)/(s_1cell.WE_DNAm.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_2cell.WE.int.reps.rem+1)/(s_2cell.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs 2cell WE', col=rgb(0,0,0,alpha=.5)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_2cell.WE.int.reps.rem+1)/(s_2cell.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_4cell.WE.int.reps.rem+1)/(s_4cell.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs 4cell WE', col=rgb(0,0,0,alpha=.5)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_4cell.WE.int.reps.rem+1)/(s_4cell.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Morula.WE.int.reps.rem+1)/(s_Morula.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs Morula WE', col=rgb(0,0,0,alpha=.5)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Morula.WE.int.reps.rem+1)/(s_Morula.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Blast.WE.int.reps.rem+1)/(s_Blast.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs Blast WE', col=rgb(0,0,0,alpha=.5)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Blast.WE.int.reps.rem+1)/(s_Blast.WE.ex+1)), pch=20, cex=.75, col='red'))
# dev.off()


# pdf(file.path(outpath, 'ExInt.Rat.1cell.rats.pdf'))
# gene.sel.MII = genes.sel.samp.rat
# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_MII.WE.int.reps.rem+1)/(s_MII.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs MII WE'))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_MII.WE.int.reps.rem+1)/(s_MII.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_1cell.WE_DNAm.int.reps.rem+1)/(s_1cell.WE_DNAm.ex+1)), pch=20, cex=.75, main='1cell WE vs 1cell DRB'))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_1cell.WE_DNAm.int.reps.rem+1)/(s_1cell.WE_DNAm.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_2cell.WE.int.reps.rem+1)/(s_2cell.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs 2cell WE'))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_2cell.WE.int.reps.rem+1)/(s_2cell.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_4cell.WE.int.reps.rem+1)/(s_4cell.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs 4cell WE'))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_4cell.WE.int.reps.rem+1)/(s_4cell.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Morula.WE.int.reps.rem+1)/(s_Morula.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs Morula WE'))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Morula.WE.int.reps.rem+1)/(s_Morula.WE.ex+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Blast.WE.int.reps.rem+1)/(s_Blast.WE.ex+1)), pch=20, cex=.75, main='1cell WE vs Blast WE'))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)/(s_1cell.WE.ex + 1)) ~ log2((s_Blast.WE.int.reps.rem+1)/(s_Blast.WE.ex+1)), pch=20, cex=.75, col='red'))
# dev.off()
# subset(gene.sel.MII, log2((s_Blast.WE.int.reps.rem+1)) > 2)

# pdf(file.path(outpath, 'Int.1cell.pdf'))
# gene.sel.MII = genes.sel.samp.rat
# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_MII.WE.int.reps.rem+1)), pch=20, cex=.75, main='1cell WE vs MII WE', col=rgb(0,0,0,alpha=.1)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_MII.WE.int.reps.rem+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_1cell.WE_DNAm.int.reps.rem+1)), pch=20, cex=.75, main='1cell WE vs 1cell DRB', col=rgb(0,0,0,alpha=.1)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_1cell.WE_DNAm.int.reps.rem+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_2cell.WE.int.reps.rem+1)), pch=20, cex=.75, main='1cell WE vs 2cell WE', col=rgb(0,0,0,alpha=.1)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_2cell.WE.int.reps.rem+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_4cell.WE.int.reps.rem+1)), pch=20, cex=.75, main='1cell WE vs 4cell WE', col=rgb(0,0,0,alpha=.1)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_4cell.WE.int.reps.rem+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_Morula.WE.int.reps.rem+1)), pch=20, cex=.75, main='1cell WE vs Morula WE', col=rgb(0,0,0,alpha=.1)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_Morula.WE.int.reps.rem+1)), pch=20, cex=.75, col='red'))

# with(genes.sel, plot(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_Blast.WE.int.reps.rem+1)), pch=20, cex=.75, main='1cell WE vs Blast WE', col=rgb(0,0,0,alpha=.1)))
# with(gene.sel.MII, points(log2((s_1cell.WE.int.reps.rem + 1)) ~ log2((s_Blast.WE.int.reps.rem+1)), pch=20, cex=.75, col='red'))
# dev.off()

