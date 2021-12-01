### INFO: Defines the clusters for the shrimp data
### DATE: 02.09.2012.
### AUTHOR:
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path='/home/members/vfranke/MyLib'
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'ScanLib.R'))

library(Cairo)
library(data.table)
#/{1} LIBRARIES


# {2} CODE

#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
rpkm.cut = 10

genome.name = 'BSgenome.Mmusculus.UCSC.mm9'

rpkm1 = 5
rpkm2 = 5
width=50

nm.lim=0
nh.lim=4
#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES


# {{3}} MAIN
genome = GenomeLoader(genome.name)

bam.files = list.files(inpath, pattern='bam$', full.names=T, recursive=T)
# bam.files = bam.files[str_detect(bam.files, 'MII') | str_detect(bam.files, 'Som')]
bam.files = bam.files[str_detect(bam.files, 's_')]

seqlen = seqlengths(genome)
seqlen = seqlen[!str_detect(names(seqlen), '_')]

stats.files = list.files(stats.path, full.names=T, pattern='RData')


l = list()
### saves the filtered bam files into a RData format
for(i in 1:length(bam.files)){
  
  l.samps = list()
  bam.file = bam.files[i]
  name = BamName(bam.file)
  print(name)
  
  cat('loading the rdata...\n')
  rdata.file = stats.files[str_detect(stats.files, name)]
  Assigner(rdata.file, 'rdata')
  cat('subsetting the rdata...\n')
  rdata.sub = subset(rdata, subset=(cut=='(20,23]' & NM <= nm.lim))
  nm = rdata.sub[,.N, by=qname]
  m = nm$N[match(rdata.sub$qname, nm$qname)]
  rdata.sub$nh = m
  rdata.sub = subset(rdata.sub, subset=nh <= 4)
  cat('getting the reads rdata...\n')
  qname.all = rdata[,length(unique(qname[NM == 0]))]
  qname = unique(rdata.sub$qname)
  
  rm(rdata, nm, m);gc()
  l.stats = list()
  library(doMC)
  registerDoMC(11)
  l.bam = foreach(chr = names(seqlen))%dopar%{
    
    print(chr)
    l.counts = list()
    which = GRanges(chr, IRanges(1, seqlen[chr]))
    param = ScanBamParam(which=which, what=c('qname'), tag='NM')
    bam = readBamGappedAlignments(bam.file, param=param)
    bam = bam[values(bam)$qname %in% qname]
    g = granges(bam)
    values(g)$nh = round(1/rdata.sub[,nh[match(values(bam)$qname, qname)]],2)
    if(length(g) != length(bam))
      stop('smth gone wrong!!!')
    g
  }
  r = unlist(GRangesList(l.bam))
  
  l = list(reads=r, total = qname.all, total.21 = length(qname))
  
  save(l, file=file.path(outpath, paste(name,'nm',nm.lim,'nh',nh.lim,'regions.RData', sep='.')))
  rm(r, l.bam,l);gc()
}
