### INFO: R Script
### DATE: 30.05.2014
### AUTHOR: Vedran Franke
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path='/home/members/vfranke/MyLib'
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

PlotClusters = function(clust.plot, name, cols='black'){
  
  if(length(cols)>1){
    annot = clust.plot$annot
    cols=cols[as.character(annot)]
  }
  
  outfile = (file.path(outpath, paste(name, 'pdf',sep='.')))
  cnts.plot = as.data.frame(values(clust.plot))[,-1]
  
  pdf(outfile, width=8, height=8)
  for(i in 1:(ncol(cnts.plot)-1)){
    for(j in (i+1):ncol(cnts.plot)){
      
      namx = colnames(cnts.plot)[i]
      namy = colnames(cnts.plot)[j]
      name = paste(namy,' VS ',namx, sep='_')
      plot(cnts.plot[[namx]], cnts.plot[[namy]], pch=20, cex=0.8, xlab=namx, ylab=namy, main=name, axes=FALSE, col=cols)
      axis(1, lwd=2)
      axis(2, lwd=2)
      if(length(cols) != 1)
        legend('bottomright', fill=cols[levels(annot)], legend=levels(annot))
    }
  }
  dev.off()
}


MergeDataTable = function(l, key=NULL, all=TRUE){
  
  if(is.null(key))
    stop('You have to give a key')
  
  l = lapply(l, function(x)setkeyv(x, cols=key))
  r = Reduce(function(x,y)merge(x,y, all=all), l)
  return(r)
}


#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES


outpath = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/RNAi/Compare_Tam-GarciaLopez'
dir.create(outpath, showWarnings=FALSE)

rdata.path = file.path(outpath, 'RData')
dir.create(rdata.path, showWarnings=FALSE)

gene.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/150523_mm9_Mus_musculus.NCBIM37.67.gtf'

reps.path = '/common/DB/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.repeatmasker.RData'

mirbase.path = '/common/DB/vfranke/Base/MirBase/mm/mm9.LO.mm10.mirbase.form.gff'

file.paths=list('/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/AccessoryData/Tam_2008_Nature/Mapped/Bowtie/n3_m20_best_strata','/common/DB/vfranke/Base/AccessoryData/GSE59254_GarciaLopez_2015_Oocyte','/common/DB/vfranke/Base/AccessoryData/GSE45983_GarciaLopez_2014_Oocyte/Mapped/Bowtie')


#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
registerDoMC(10)

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

mirbase = read.table(mirbase.path, header=FALSE, sep='\t')
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
bamfiles = bamfiles[c(1,2,4,5)]

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
  seqlevels(reads, force=TRUE) = setdiff(seqlevels(reads), 'chrM')
  lg[[name]]=reads
}

lbam = list()
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
names(lbam)[3:4] = paste(c('GL15','GL14'), names(lbam)[3:4], sep='.')
names(lbam) = str_replace(names(lbam),'.cl.*','')
names(lbam)[1:2] = paste('Tam', str_replace(names(lbam)[1:2],'oocytesmallRNA_',''), sep='.')

stats = do.call(rbind, lapply(lbam, '[[', 'stats'))

rownames(stats) = names(lbam)
write.table(stats, file.path(outpath, paste(Dater(),'OOSmallStats.txt',sep='')), row.names=T, col.names=T, quote=F, sep='\t')

sannot = do.call(rbind, lapply(lbam, '[[', 'annot'))
sannot = cbind(data.table(id=rep(names(lbam), each=2)), sannot)
write.table(sannot, file.path(outpath, paste(Dater(),'OOSmallAnnotStats.txt',sep='')), row.names=F, col.names=T, quote=F, sep='\t')

r = reduce(unlist(GRangesList(lapply(lbam, '[[', 'clust'))))
rclust = suppressWarnings(reduce(resize(r, width=width(r)+width, fix='center')))
rclust = resize(rclust, width= width(rclust)-width, fix='center')

lr = list()
lr = foreach(i = 1:length(lbam))%dopar%{
  
  print(i)
  name = names(lbam)[i]
  set = lbam[[i]]
  clust = set$clust
  fo = data.table(as.matrix(findOverlaps(clust,rclust)))
  fo$Uniq = clust$Uniq[fo$queryHits]
  fo$Mult = clust$Mult[fo$queryHits]
  fo = fo[,lapply(.SD, sum), by='subjectHits']
  fo$queryHits=NULL
  setnames(fo, 2:3, paste(name, names(fo)[2:3],sep='.'))
  return(fo)
}
m = MergeDataTable(lr, 'subjectHits')
m[is.na(m)] = 0
m = m[,!str_detect(colnames(m),'19to24'),with=FALSE]
m.uniq = m[,str_detect(colnames(m),'Uniq'), with=FALSE]
m.mult = m[,str_detect(colnames(m),'Mult'), with=FALSE]
tot.cnts =  sapply(lg, function(x)length(unique(names(x))))
tot.cnts = tot.cnts[!str_detect(names(tot.cnts),'19to24')]

m.uniq.rpkm = rpkm(m.uniq, gene.length=width(rclust), lib.size = tot.cnts)

# rpkm.border = quantile(Reduce(pmax, as.data.frame(m.uniq.rpkm)), seq(0,1,.01))['95%']
# clust.ind = rowSums(m.uniq.rpkm >= rpkm.border) > 0
# clust.ind = Reduce(pmax, as.data.frame(m.uniq.rpkm))
clust.ind = rowMeans(as.data.frame(m.uniq.rpkm))
clust.ind = head(order(-clust.ind), 2000)
clust.sel = rclust[clust.ind]
m.cnts = m.uniq + m.mult
m.cnts.sel = m.cnts[clust.ind,]
m.rpkm = as.data.frame(log10(rpkm(as.data.frame(m.cnts.sel+1), gene.length=width(clust.sel), lib.size=tot.cnts)))
m.rpkm[m.cnts.sel == 0] = 0
m.cpm = as.data.frame(log10(t(t(as.data.frame(m.cnts.sel+1))*(1e7/tot.cnts))))
m.cpm[m.cnts.sel == 0] = 0


clust.sel$annot = AnnotateRanges(clust.sel, lannot)
values(clust.sel) = cbind(values(clust.sel), as.data.frame(m.cpm))
clust.sel = clust.sel[!clust.sel$annot %in% c('tRNA','rRNA.ens','rRNA.repmask')]
clust.sel$annot[clust.sel$annot %in% c('LINE','SINE','LTR')] = 'Repeats'
clust.sel$annot[clust.sel$annot == 'linc.exon'] = 'exon'
clust.sel$annot[clust.sel$annot == 'linc.exon.anti'] = 'exon.anti'
clust.sel$annot = factor(clust.sel$annot)

cols = c('chartreuse3','darkorange','cornflowerblue','darkorchid3','black','gray','firebrick2')
names(cols) = levels(clust.sel$annot)


PlotClusters(clust.sel, paste(Dater(),'OOsmallRNA.ScatterPlot', length(clust.sel), sep='.'), cols)
PlotClusters(clust.sel, paste(Dater(),'OOsmallRNA.ScatterPlot.black', length(clust.sel), sep='.'), 'black')

clust.sel.woe = clust.sel[!(str_detect(clust.sel$annot, 'exon') | clust.sel$annot == 'MT')]
clust.sel.woe$annot = factor(as.character(clust.sel.woe$annot))
PlotClusters(clust.sel.woe, paste('OOsmallRNA.ScatterPlot.WoExons_MT', rpkm.border, sep='.'), cols)

clust.sel.woe = clust.sel[!(str_detect(clust.sel$annot, 'exon') | clust.sel$annot == 'Repeats')]
clust.sel.woe$annot = factor(as.character(clust.sel.woe$annot))
PlotClusters(clust.sel.woe, paste('OOsmallRNA.ScatterPlot.WoExons_Repeats', rpkm.border, sep='.'), cols)

# correlations
cor(as.data.frame(values(clust.sel))[,-1], method='spearman')



#/{{3}} MAIN

#/{2} CODE
