### INFO: Mapping statistics on the Dicer data
### DATE: 05.07.2012.
rm(list=ls());gc()
### AUTHOR: Vedran Franke


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path='/home/members/vfranke/MyLib'
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'ChipSeq.Functions.R'))
library(Rsamtools)
library(doMC)
library(data.table)
#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS
#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES
inpath = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/MappingStats/ShrimpStrata/RDataStats'

outpath='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/SuppFigure10/'
dir.create(outpath, showWarnings=F)

sampnames.path='/home/members/vfranke/Projects/Dicer_RNASeq_29062012/Scripts/SampleLabels.R'

bam.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/Mapped/ShrimpStrata'

stats.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/MappingStats/ShrimpStrata/RData4'

regions.path='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData/SingleRegions/regions.K.txt'

mosir.path='/common/USERS/vfranke/Base/Genomes/mm9/Sequence/Chrs/pCAG-EGFP_MosIR.fasta'
#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
registerDoMC(6)

s = 'TTGTGGCTGCCAGCACGCGCACGCCCGAAGACTCCAACAGCCTAGGTACCATAA'
# registerDoMC(9)
#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES
source(sampnames.path)
mosir.seq = DNAString(as.character(read.DNAStringSet(mosir.path)))
m1 = matchPattern(s, mosir.seq)
m2 = matchPattern(reverseComplement(DNAString(s)), mosir.seq)

bam.files = list.files(bam.path, pattern='bam$', full.names=T, recursive=T)
bam.files = bam.files[str_detect(bam.files, 's_')]

mosir = GRanges('pCAG-EGFP_MosIR', IRanges(1,6660))


registerDoMC(9)
l.bam = foreach(i = 1:length(bam.files))%dopar%{
  
  l.samps = list()
  bam.file = bam.files[i]
  name = BamName(bam.file)
  print(name)
  
  rdata.file = list.files(stats.path, pattern=name, full.names=T)
  Assigner(rdata.file, 'rdata')
  qname = rdata[,qname[cut== '(20,23]' & NM <= 0]]
  
  param = ScanBamParam(which=mosir, what=c('qname'))
  bam = readBamGappedAlignments(bam.file, param=param)
  bam = bam[values(bam)$qname %in% qname]
  list(bam=bam, total = rdata[,length(unique(qname[NM == 0]))], total.23 = length(unique(qname)))
}
names(l.bam) = BamName(bam.files)
save(l.bam, file=file.path(outpath, 'MosReads.RData'))
bams = lapply(l.bam, '[[', 'bam')
# bams = list(T3T_DcrO = c(bams[[1]], bams[[2]], bams[[3]]),
# T3T_DcrS = c(bams[[4]], bams[[5]], bams[[6]]),
# T3T_Utra = c(bams[[7]], bams[[8]], bams[[9]]))
bams = lapply(bams, function(x){
  g = granges(x)
  sind = as.character(strand(g)) == '+'
  list(covp = coverage(g[sind])[[23]], covm = coverage(g[!sind])[[23]])})
total = sapply(l.bam , '[[', 'total')
# total = c(sum(total[1:3]), sum(total[4:6]), sum(total[7:9]))
total = 1e6/total

for(i in 1:length(bams)){
  
  name = names(bams)[i]
  print(name)
  cov = bams[[i]]
  cov = lapply(cov, function(x)as.vector(x*total[i]))
  cov[[2]] = -cov[[2]]
  ylim = trunc(range(sapply(cov, range)))+c(-1,1)
  xlim = c(2530,3010)
  lwd = 4
  CairoPNG(file.path(outpath, paste(name, 'red.png', sep='')), width=800, height=600)
  plot(cov$covp, type='l', xlim=xlim, ylim=ylim,axes=F)
  polygon(cov$covp, col='black')
  lines(cov$covm)
  polygon(cov$covm, col='black')
  if(name == "s_ES_WT_tran"){
    st = start(m2):end(m2)
    lines(st, rep(0, length(st)), col='red', lwd=5)
  }
  axis(2, at=ylim, lwd=lwd)
  # axis(1, at = seq(xlim[1], xlim[1]+ 50))
  dev.off()
  
}
