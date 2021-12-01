#!/common/WORK/vfranke/bin/R/R-3.1.0/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path='/home/members/vfranke/MyLib'
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'ScanLib.R'))
library(data.table)
library(stringr)
library(doMC)
library(GenomicAlignments)
library(genomation)

#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS

#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES
file = '/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_1cell.WE/s_1cell.WE.bam'
file = %FILE


outpath = %OUTPATH

r.outpath = file.path(outpath,'Count_ReadsInExons_RData');
dir.create(r.outpath, showWarnings=FALSE)

mm9.gtf.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/22012012_EnsemblGenes_BioMart.form.gtf'

ens.genes.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/25012012.EnsemblGenesAll_BioMart.form.txt'

rep.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.repeatmasker.RData'

encode.path = '/common/DB/vfranke/Base/Encode/mm9/RNASeq/Mapped'


ooreg.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/GenomeActivation/OOcyteRegions/OOcyteRegions.lower.5.bed'
#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
registerDoMC(21)
#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES

# {{3}} MAIN

# reads in the ensembl annotation
gtf = read.table(mm9.gtf.path, header=FALSE, sep='\t')
gtf = gtf[!str_detect(gtf[,1],'NT'),]
gtf[gtf[,1] == 'chrMT',1] = 'chrM'

gtf = GffToGRanges(gtf, 'exon')
gids = data.table(as.data.frame(unique(values(gtf)[c('gene_id','transcript_id','exon_id')])))
annot = read.table(ens.genes.path, header=TRUE, sep='\t')
gids$biotype = annot$biotype[match(gids$gene_id, annot$ens.gene.id)]
gids = gids[gids$biotype %in% c('protein_coding','lincRNA'),]
gids = gids[!duplicated(gids$exon_id)]

exon = gtf[gtf$exon_id %in% gids$exon_id]
values(exon) = values(exon)$exon_id
exon = exon[!duplicated(exon$value)]

exon.unique = exon[countOverlaps(exon, exon, ignore.strand=TRUE) == 1]

name = BamName(file)
print(name)
chrs = chrFinder(file)
chrs = chrs[chrs$chr != 'chrM',]

cnts = foreach(chr = chrs$chr)%dopar%{
  # cnts = foreach(chr = chrs$chr[1:2])%dopar%{
  
  
  print(chr)
  # reads in the reads
  cat('Reading files...\n')
  which.ranges = GRanges(chr, IRanges(1, chrs$chr.len[chrs$chr == chr]))
  if(str_detect(name, 'SE')){
    bam = readGAlignmentsFromBam(file, param=ScanBamParam(which=which.ranges, tag='NH'), use.names=TRUE)
    bams = bam
  }else{
    bam = suppressWarnings(readGAlignmentPairs(file, param=ScanBamParam(which=which.ranges), use.names=TRUE))
    bams = readGAlignmentsFromBam(file, param=ScanBamParam(which=which.ranges, tag='NH'), use.names=TRUE)
  }
  
  gbam = grglist(bam)
  tab = unique(data.table(names=names(bams), V1=values(bams)$NH))
  tab = tab[tab$names %in% names(bam)]
  rm(bam, bams);gc()
  
  # counts the reads in exons
  cat('Counting reads in models...\n')
  foe = data.table(as.matrix(findOverlaps(gbam, exon.unique, ignore.strand=TRUE)))
  foe$exon_id = exon.unique$value[foe$subjectHits]
  
  foe.n = foe[, length(unique(queryHits)), by=exon_id]
  setnames(foe.n,2,name)
  
  return(foe.n)
}
cnts = cnts[sapply(cnts, function(x)any(class(x)  %in% 'data.table'))]
tab = rbindlist(cnts)


cat('Saving data...\n')
save(tab, file = file.path(r.outpath, paste(name,'UniqueExon_Counts.RData',sep='.')))

cat('Thank you, and goodbye!...\n')
