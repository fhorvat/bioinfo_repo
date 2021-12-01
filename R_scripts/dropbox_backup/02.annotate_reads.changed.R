### INFO: R Script
### DATE: 30.05.2014
### AUTHOR: Vedran Franke
rm(list=ls());gc()
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/Analysis/2020_paper/small_RNA_reads_annotation/vfranke_test")
options(width = 185)

# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path='/home/members/vfranke/MyLib/RFun'
source(file.path(lib.path, 'ScanLib.R'))
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'FormatConverters.R'))

library(data.table)
library(stringr)
library(ggplot2)
# library(Cairo)
library(doMC)
library(GenomicAlignments)
library(reshape2)
library(edgeR)


#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS

#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES


outpath = getwd()

rdata.path = file.path(outpath, 'RData')
dir.create(rdata.path, showWarnings=FALSE)

gene.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/150523_mm9_Mus_musculus.NCBIM37.67.gtf'

reps.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.RepeatMasker.RData'

mirbase.path = '/common/DB/vfranke/Base/MirBase/mm/mm9.LO.mm10.mirbase.form.gff'

file.paths=list('/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/AccessoryData/Tam_2008_Nature/Mapped/Bowtie/n3_m20_best_strata',
                '/common/DB/vfranke/Base/AccessoryData/GSE59254_GarciaLopez_2015_Oocyte',
                '/common/DB/vfranke/Base/AccessoryData/GSE45983_GarciaLopez_2014_Oocyte/Mapped/Bowtie')


#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
# registerDoMC(10)

width=100

genome.name = 'BSgenome.Mmusculus.UCSC.mm9'
#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES


# {{3}} MAIN

Assigner(reps.path, 'reps')
genes = fread(gene.path)
genes = genes[!str_detect(genes$V1,'NT')]
genes = genes[!str_detect(genes$V1,'MT')]
genes$V1 = paste('chr', genes$V1, sep='')
genes$V9 = str_replace_all(genes$V9,'"','')
gff = GffToGRanges(as.data.frame(genes), 'exon')
exon = gff[gff$gene_biotype == 'protein_coding']
intron = unlist(range(split(exon, exon$transcript_id)))
pseudo = gff[gff$gene_biotype == 'pseudogene']
exon.anti = exon
strand(exon.anti) = ifelse(strand(exon.anti) == '+','-','+')

linc.exon = gff[gff$gene_biotype == 'lincRNA']
linc.intron = unlist(range(split(linc.exon, exon$transcript_id)))
linc.exon.anti = linc.exon
strand(linc.exon.anti) = ifelse(strand(linc.exon.anti) == '+','-','+')

rRNA.ens = gff[gff$gene_biotype == 'rRNA']
rRNA.rep = reps[reps$repClass == 'rRNA']

tRNA = reps[reps$repClass == 'tRNA']

mirbase = read.table(mirbase.path, header=FALSE, sep='\t', stringsAsFactors = F)
mirbase = na.omit(mirbase)
mirbase = GffToGRanges(mirbase)

lreps = list(MT = reps[str_detect(reps$repName, 'MT[ABCD]')],
             SINE = reps[reps$repClass == 'SINE'],
             LINE = reps[reps$repClass == 'LINE'],
             LTR = reps[reps$repClass == 'LTR'])
lannot1 = list(rRNA.ens = rRNA.ens,
               rRNA.repmask = rRNA.rep,
               tRNA = tRNA,
               mirna                =       mirbase,
               pseudo               =       pseudo)
lannot2 = list(exon.anti        =       exon.anti,
               linc.exon.anti = linc.exon.anti,
               exon         =       exon,
               linc.exon    =       linc.exon)
lannot = c(lannot1, lreps, lannot2)
lannot = GRangesList(lapply(lannot, function(x){values(x)=NULL;x}))

genome = GenomeLoader(genome.name)

bamfiles = unlist(lapply(file.paths, function(x)list.files(x, full.names=TRUE, recursive=TRUE, pattern='bam$')))
bamfiles = bamfiles[!str_detect(bamfiles, 'Bowtie2')]
bamfiles = bamfiles[c(1,2)]

lg = list()
for(i in 1:length(bamfiles)){
  
  bamfile = bamfiles[i]
  name = BamName(bamfile)
  name = str_replace(name, '_genome', '')
  name = str_replace(name, '.cleaned.fq', '')
  if(str_detect(bamfile,'GarciaLopez_2015'))
    name = 'GL15.Oocyte'
  if(str_detect(bamfile,'GarciaLopez_2014'))
    name = 'GL14.Oocyte'
  print(name)
  
  reads = readGAlignments(bamfile, use.names=TRUE, param=ScanBamParam(tag='NM'))
  seqlevels(reads, pruning.mode="coarse") = setdiff(seqlevels(reads), 'chrM')
  lg[[name]]=reads
}

lbam = list()
for(i in 1:length(bamfiles)){
  
  i <- 2
  
  bamfile = bamfiles[i]
  name = BamName(bamfile)
  name = str_replace(name, '_genome', '')
  name = str_replace(name, '.cleaned.fq', '')
  if(str_detect(bamfile,'GarciaLopez_2015'))
    name = 'GL15.Oocyte'
  if(str_detect(bamfile,'GarciaLopez_2014'))
    name = 'GL14.Oocyte'
  print(name)
  
  reads = lg[[i]]
  g = granges(reads, use.mcols=FALSE)
  g$name = paste(name, names(g), sep='.')
  
  gannot = AnnotateRanges(g, lannot)
  tab = data.table(name=names(reads))
  tab = tab[,.N, by=name]
  tab$width = width(reads)[match(tab$name, names(reads))]
  tab$short = tab$width >= 21 & tab$width <= 23
  tab$uniq = ifelse(tab$N==1, 'Uniq', 'Mult')
  tab$NM = values(reads)$NM[match(tab$name, names(reads))]
  tab = tab[tab$N <= 5]
  
  tab$annot = gannot[match(tab$name, names(reads))]
  tab$annot = factor(tab$annot, levels=c(names(lannot),'None'), ordered=TRUE)
  tabm = tab[order(tab$annot),]
  tabm = tabm[!duplicated(tabm$name)]
  dtabm = tabm[tabm$short,length(unique(name)),by=list(annot, uniq)]
  dtabm = dcast.data.table(uniq~annot, data=dtabm)
  
  g.short = g[names(g) %in% tab$name[tab$short & tab$NM <=1]]
  clust = reduce(resize(g.short, width=width(g.short)+width, fix='center'))
  clust = resize(clust, width= width(clust)-width, fix='center')
  
  fo = data.table(as.matrix(findOverlaps(clust, g.short)))
  fo$rname = names(g.short)[fo$subjectHits]
  fo$uniq   = tab$uniq[match(fo$rname, tab$name)]
  
  foc = fo[,length(unique(subjectHits)),by=c('queryHits','uniq')]
  fos = dcast.data.table(foc, queryHits~uniq,  value.var='V1', fun=sum)
  fos[is.na(fos)] = 0
  clust$Uniq = 0
  clust$Uniq[fos$queryHits] = fos$Uniq
  clust$Mult = 0
  clust$Mult[fos$queryHits] = fos$Mult
  setnames(fos, 2:3, paste(name,colnames(fos)[2:3],sep='.'))
  
  stats = tab[,length(unique(name)), by=c('short','uniq') ]
  stats = stats[,list(tot=sum(V1), tot.uniq=sum(V1[uniq=='Uniq']), tot.21=sum(V1[short]),tot.21.uniq = sum(V1[short& uniq=='Uniq']))]
  rownames(stats) = name
  l=list(clust=clust, stats=stats, annot=dtabm)
  save(l, file = file.path(rdata.path, paste(name, 'RData', sep='.')))
  lbam[[name]] = l
}
names(lbam) = str_replace(BamName(bamfiles),'_genome','')
# names(lbam)[3:4] = paste(c('GL15','GL14'), names(lbam)[3:4], sep='.')
names(lbam) = str_replace(names(lbam),'.cl.*','')
names(lbam)[1:2] = paste('Tam', str_replace(names(lbam)[1:2],'oocytesmallRNA_',''), sep='.')

stats = do.call(rbind, lapply(lbam, '[[', 'stats'))

rownames(stats) = names(lbam)
write.table(stats, file.path(outpath, paste(Dater(),'OOSmallStats.txt',sep='')), row.names=T, col.names=T, quote=F, sep='\t')

sannot = lapply(lbam, '[[', 'annot') %>% dplyr::bind_rows(.)
sannot = cbind(data.table(id=rep(names(lbam), each=2)), sannot)
write.table(sannot, file.path(outpath, paste(Dater(),'OOSmallAnnotStats.txt',sep='')), row.names=F, col.names=T, quote=F, sep='\t')

