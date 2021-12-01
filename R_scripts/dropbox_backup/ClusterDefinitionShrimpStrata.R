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
# {{1}} FUNCTIONS
# -------------------------------------- #
#
GetOverlaps = function(reg1,reg2, colname=NULL){
  
  if(is.null(colname))
    stop('Colname needs to be defined')
  fo = data.table(as.matrix(findOverlaps(reg1,reg2)), ignore.strand=T)
  fo$weight = values(reg2)[[colname]][fo$subjectHits]
  fo = fo[,sum(weight), by=queryHits]
  v = rep(0, length(reg1))
  v[fo$queryHits] = fo$V1
  return(v)
}


# -------------------------------------- #
# finds the clusters in short RNA data
FindClusters = function(reads, total, rpkm1 = 5, rpkm2 = 10, clust.width=100, e=1e6){
  
  cat('Reducing the regions...\n')
  dregs = reduce(reads, ignore.strand=T)
  
  cat('clustnum:', length(dregs), '\n')
  co = GetOverlaps(dregs, reads,'nh')
  cat('First step filtering...\n')
  rpkm = (co )*(e/total)
  rind = rpkm > rpkm1
  
  cat('Distance merging...\n')
  wregs = dregs[rind]
  # wregs = dregs
  values(wregs)$counts = co[rind]
  values(wregs)$rpkm = rpkm[rind]
  # eregs = suppressWarnings(resize(wregs, width= width(wregs)+clust.width, fix='center'))
  # fo = as.matrix(findOverlaps(eregs, eregs, ignore.strand = T))
  # dup.ind = which(fo[-nrow(fo),1] == fo[-1,2] & fo[-1,1] == fo[-(nrow(fo)),2] )
  # fo = fo[-dup.ind, ]
  
  ### cluster graphs
  cat('Graph clustering...\n')
  # m = ClusterGraphs(fo)
  a.regs = reduce(resize(wregs, width=width(wregs) + clust.width, fix='center'), ignore.strand=T)
  a.regs = resize(a.regs, width=width(a.regs) - clust.width, fix='center')
  values(a.regs)$counts = GetOverlaps(a.regs, reads,'nh')
  
  
  cat('Cluster summarization...\n')
  # m = m[order(m[,1]),]
  # values(wregs)$cwidth = width(wregs)
  # cregs = split(wregs, m[,2])
  # el = elementLengths(cregs) == 1
  # s.regs = cregs[el]
  # w.regs = GRangesList(lapply(cregs[!el], function(x){
  # g = GRanges(unique(as.character(seqnames(x))), IRanges(min(start(x)), max(end(x))))
  # values(g)$counts = sum(values(x)$counts)
  # values(g)$rpkm = sum(values(x)$rpkm)
  # values(g)$cwidth = sum(width(x))
  # g}))
  
  
  # a.regs = c(s.regs, w.regs)
  # a.regs = a.regs[order(as.numeric(names(a.regs)))]
  # a.regs = unlist(a.regs)
  # values(a.regs)$rpkm.c = (values(a.regs)$counts / values(a.regs)$cwidth)*(e/total)
  values(a.regs)$rpkm.c = (values(a.regs)$counts)*(e/total)
  
  cat('Second step filtering...\n')
  # f.ind = values(a.regs)$rpkm.c > rpkm2
  regs.sel = a.regs
  # cregs.sel = cregs[f.ind]
  
  cat('Returning the clusters...\n')
  return(reg.sel = regs.sel)
  # return(list(regs.sel = regs.sel, reg.sel.list = cregs.sel))
}


# -------------------------------------- #
# designates the elements in the graph based on the connectivity matrix
ClusterGraphs = function(fo){
  
  n = unique(fo[,1])
  v = vector()
  lc = list()
  for(i in n){
    
    cat(i,'\r')
    if(i %in% v)
      next
    
    p = fo[fo[,1] == i,2]
    s = vector()
    while(length(p) > 0){
      
      for(j in p){
        s = c(s, j)
        p = unique(c(p, fo[fo[,1] == j,2], fo[fo[,2] == j,1]))
      }
      p = p[!p %in% s]
    }
    # s = unique(c(i,s))
    v = unique(c(v, s))
    lc = c(lc, list(s))
  }
  m = cbind(unlist(lc), rep(1:length(lc), times=sapply(lc, length)))
  return(m)
}

# m1 = matrix(c(1:5,1,3,4,4,5), ncol=2)
# m2 = matrix(c(1:5,1,1,3,5,5), ncol=2)
# -------------------------------------- #
# merges sets of cluster
ClusterMerge = function(regs1, regs2, method='intersect', ignore.strand=T){
  
  cat('method:', method,'\n')
  if(method == 'intersect'){
    fo = as.matrix(findOverlaps(regs1, regs2))
    o = reduce(c(regs1[fo[,1]], regs2[fo[,2]]), ignore.strand=ignore.strand)
  }
  if(method == 'union'){
    o = reduce(c(regs1, regs2), ignore.strand=ignore.strand)
    
  }
  
  return(regs.sel=o)
}


# -------------------------------------- #
# get annotation
GetAnnotation = function(l.paths){
  
  mirna.path = l.paths$mirna.path
  exons.path = l.paths$exons.path
  genes.path = l.paths$ens.genes.path
  reps.path = l.paths$repeat.path
  
  mirna = read.table(mirna.path, header=F)
  mirna[,9] = str_replace_all(mirna[,9],'=',' ')
  mirna = mirna[!(is.na(mirna[,4]) | is.na(mirna[,5])),]
  mirna = GffToGRanges(mirna)
  
  exons = read.table(exons.path, header=F, sep='\t')
  exons = GffToGRanges(exons)
  
  genes = read.table(genes.path, header=T)
  genes = BedToGRanges(genes, values=T)
  genes = genes[!values(genes)$biotype %in% c('pseudogene','miRNA')]
  exons = exons[values(exons)$gene_id %in% values(genes)$ens.gene.id]
  
  gr = reduce(genes)
  ex = reduce(exons)
  
  Assigner(reps.path, 'reps')
  reps = updateObject(reps)
  reps = reps[!str_detect(as.vector(seqnames(reps)), 'random')]
  seqlevels(reps) = unique(as.character(seqnames(reps)))
  reps = reduce(reps)
  
  list(
    miRNA=mirna,
    Genes=ex,
    Repeats=reps)
  
}


# -------------------------------------- #
# annotates each element of the regions with the corresponding overlapping annotations
AnnotateRegions = function(regs, l.annot){
  
  d = do.call(cbind, lapply(l.annot, function(x)suppressWarnings(countOverlaps(regs, x))))
  d[d>1] = 1
  d = t(t(d) * 1:ncol(d))
  d[d == 0] = max(d)+1
  da = t(data.frame(d))
  mi = do.call(pmin, split(da, 1:ncol(d)))
  ind = c(names(l.annot), 'Other')[mi]
  ind
}

# -------------------------------------- #
PlotClusters = function(expr1, expr2,ind , outpath, name, xlab='', ylab='',mcw=NULL, cols=NULL){
  
  lwd=5
  cex.axis=2
  cex=1.5
  ind = as.factor(ind)
  if(is.null(cols)){
    cols = c('red','cornflowerblue','black','darkgray')
  }
  outname.sq = paste(name,'mcw',mcw, 'sq.png', sep='.')
  CairoPNG(file.path(outpath, outname.sq), width=800, height=800)
  plot(expr1, expr2, pch=20, xlab=xlab, ylab=ylab, col = cols[as.numeric(ind)], axes=F, xlim=c(0,6), ylim=c(0,6), cex=cex)
  axis(1, at=0:5, lwd=lwd, cex=cex.axis, tcl=-1, labels=F,  padj=5)
  axis(2, at=0:5, lwd=lwd, cex=cex.axis, tcl=-1, labels=F,  padj=5)
  legend('topright', legend=levels(ind), fill=cols, bty='n')
  dev.off()
  
  outname.re = paste(name,'mcw',mcw, 're.png', sep='.')
  CairoPNG(file.path(outpath, outname.re), width=1000, height=600)
  plot(expr1, expr2, pch=20, xlab=xlab, ylab=ylab, col = cols[as.numeric(ind)], axes=F, xlim=c(0,6), ylim=c(0,6), cex=cex)
  axis(1, at=0:5, lwd=lwd, cex=cex.axis, tcl=-1, labels=F)
  axis(2, at=0:5, lwd=lwd, cex=cex.axis, tcl=-1, labels=F)
  legend('topright', legend=levels(ind), fill=cols, bty='n')
  dev.off()
}

#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES
inpath='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/ClusterDefinitionShrimp'

outdir = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/ClusterDefinitionShrimp'
outpath = file.path(outdir, OutpathDate('Clust'))
dir.create(outpath, showWarnings=F)

plot.path = file.path(outpath, 'Plots')
dir.create(plot.path, showWarnings=F)

clust.path = file.path(outpath, 'Clust')
dir.create(clust.path, showWarnings=F)

stats.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/MappingStats/ShrimpStrata/RData4'


acc.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData'

watanabe.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData/Watanabe_2008_Nature/Watanabe_endoSirna.mm9.bed'

l.paths = list(
  mirna.path = '/common/USERS/vfranke/Base/MirBase/mm9.LO.mm10.mirbase.form.gff',
  ens.genes.path = '/common/USERS/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/25012012.EnsemblGenesAll_BioMart.form.txt',
  repeat.path='/common/USERS/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.repeatmasker.RData',
  exons.path = '/common/USERS/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/22012012_EnsemblGenes_BioMart.form.gtf'
)

regions.path='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData/SingleRegions/regions.K.txt'


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


### reads in the count matrices
dat = gsub('\\s+','-',date())
dat = gsub(':','_',dat)

l.files = list.files(inpath, full.names=T, pattern='regions')
l.files = l.files[str_detect(l.files, 's_')]
l.files = l.files[str_detect(l.files, 'nh')]
l.mat = list()
for(i in 1:length(l.files)){
  
  l.file = l.files[i]
  name = str_replace(basename(l.file),'.nm.+$','')
  print(name)
  Assigner(l.file, 'l')
  l.mat[[name]] = l
}

bams = lapply(l.mat, function(x)x$reads)
rpkm1=3
rpkm2=3
width=100
s4.clust = FindClusters(bams[["s_ES_Som.Dicer_Rsc_utra_r2"]], l.mat[["s_ES_Som.Dicer_Rsc_utra_r2"]]$total, rpkm1, rpkm2, width,1e6)
s11.clust = FindClusters(bams[["s_ES_Som.Dicer_Rsc_utra_r1"]], l.mat[["s_ES_Som.Dicer_Rsc_utra_r1"]]$total, rpkm1, rpkm2, width,1e6)

rpkm1=3
rpkm2=3
width=50
o2.clust = FindClusters(bams[["s_ES_MII.Dicer_Rsc_utra_r1"]], l.mat[["s_ES_MII.Dicer_Rsc_utra_r1"]]$total, rpkm1, rpkm2, width,1e6)
o3.clust = FindClusters(bams[["s_ES_MII.Dicer_Rsc_utra_r2"]], l.mat[["s_ES_MII.Dicer_Rsc_utra_r2"]]$total, rpkm1, rpkm2, width,1e6)

length(s4.clust)
length(s11.clust)
length(o2.clust)
length(o3.clust)

o.clust = ClusterMerge(o2.clust, o3.clust)
s.clust = ClusterMerge(s4.clust, s11.clust)

length(o.clust)
length(s.clust)

d.cnts = do.call(cbind, lapply(bams, function(x)GetOverlaps(o.clust, x, 'nh')))
total = unlist(lapply(l.mat, '[[', 'total'))
d.cnts.norm = data.frame(t(t(d.cnts)/total*1e6))

regions = BedToGRanges(read.table(regions.path, header=T), values=T)
regions = regions[!str_detect(values(regions)$name, 'Custom')]
countOverlaps(regions, o.clust)
# GetOverlaps(regions[2], bams[["s_ES_MII.Dicer_Rsc_utra_r1"]], 'nh')
# GetOverlaps(regions[2], bams[["s_ES_MII.Dicer_Rsc_utra_r2"]], 'nh')

# opt = o.clust[countOverlaps(o.clust,regions[2]) == 1]
# g1=GetOverlaps(opt, bams[["s_ES_MII.Dicer_Rsc_utra_r1"]], 'nh')
# g1*1e6/l.mat[["s_ES_Som.Dicer_Rsc_utra_r1"]]$total
# g2=GetOverlaps(opt, bams[["s_ES_MII.Dicer_Rsc_utra_r2"]], 'nh')
# g2*1e6/l.mat[["s_ES_Som.Dicer_Rsc_utra_r2"]]$total
# ------------------------------------------------------------------- #
# 1. O table
cp = GetOverlaps(o.clust,bams[['s_ES_MII.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_MII.Dicer_Rsc_utra_r2']]) == '+'])
cm = GetOverlaps(o.clust,bams[['s_ES_MII.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_MII.Dicer_Rsc_utra_r2']]) == '-'])
perc = cp/(cp+cm)

do = data.frame(as.data.frame(o.clust), d.cnts.norm)
do$strand = '*'
do$strand[perc >.8] = '+'
do$strand[perc <.2] = '-'

l.annot = GetAnnotation(l.paths)
do$class = AnnotateRegions(o.clust, l.annot)
write.table(do, file.path(clust.path, paste('O.combined',dat,'txt',sep='.')), row.names=F, col.names=T, quote=F, sep='\t')


# 1.2 mirnas
mirna = read.table(mirna.path, header=F)
mirna[,9] = str_replace_all(mirna[,9],'=',' ')
mirna = mirna[!(is.na(mirna[,4]) | is.na(mirna[,5])),]
mirna = GffToGRanges(mirna)

mov = data.table(as.matrix(findOverlaps(o.clust, mirna)))
mo = mov
mo$id = values(mirna)$X_ID[mov$subjectHits]
id = mo[,paste(unique(id), collapse=':'), by=queryHits]
dom = cbind(do[unique(id$queryHits),], mirna=id$V1)
dom = dom[order(-dom$s_ES_MII.Dicer_Rsc_utra_r2),]
write.table(head(dom,60), file.path(clust.path, 'mirna.60.txt'), row.names=F, col.names=T, quote=F, sep='\t')


# 1.3 S table

# ------------------------------------------------------------------- #
# 2. O - S Cluster Plotting
all.clust = ClusterMerge(o.clust, s.clust, method='union')
d.cnts.all = do.call(cbind, lapply(bams, function(x)GetOverlaps(all.clust, x, 'nh')))
total = unlist(lapply(l.mat, '[[', 'total'))


d.all.norm = data.frame(t(t(d.cnts.all)/total*1e6))
annot.all = AnnotateRegions(all.clust, l.annot)
d.all = as.data.frame(all.clust)
d.all$class = annot.all


cp.all.o = GetOverlaps(all.clust,bams[['s_ES_MII.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_MII.Dicer_Rsc_utra_r2']]) == '+'], 'nh')
cm.all.o = GetOverlaps(all.clust,bams[['s_ES_MII.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_MII.Dicer_Rsc_utra_r2']]) == '-'], 'nh')
perc.o = cp.all.o/(cp.all.o+cm.all.o)
strand.o = rep('*', nrow(d.all))
strand.o[perc.o > 0.8] = '+'
strand.o[perc.o < 0.2] = '-'

cp.all.s = GetOverlaps(all.clust,bams[['s_ES_Som.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_Som.Dicer_Rsc_utra_r2']]) == '+'], 'nh')
cm.all.s = GetOverlaps(all.clust,bams[['s_ES_Som.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_Som.Dicer_Rsc_utra_r2']]) == '-'], 'nh')
perc.s = cp.all.s/(cp.all.s+cm.all.s)
strand.s = rep('*', nrow(d.all))
strand.s[perc.s > 0.8] = '+'
strand.s[perc.s < 0.2] = '-'


d.all = data.frame(d.all[,-5], strand.O=strand.o, strand.S = strand.s, d.all.norm)

write.table(d.all, file.path(clust.path, paste('O - S.combined.expr',dat,'txt',sep='.')), row.names=F, col.names=T, quote=F, sep='\t')
# d.all = read.table(file.path(clust.path, 'O - S.combined.expr.txt'),header=T)

# ------------------------------------------------------------------- #
# with mirna
# plot.path = file.path(outpath, 'PlotsTest')
dir.create(plot.path, showWarnings=F)
#### utra all
# d.all = read.table('/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/ClusterDefinitionShrimp/2012-10-02_Clust/Clust/O.combined.Swo.txt', header=T)
PlotClusters(expr1=log10(d.all$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.all$s_ES_Som.Dicer_Rsc_utra_r2+1),
             ind = d.all$class,
             outpath=plot.path,
             name='O3.vs.S11.utra',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.all$width))



PlotClusters(expr1=log10(d.all$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.all$s_ES_WT_utra+1),
             ind = d.all$class,
             outpath=plot.path,
             name='O3.vs.WT.utra',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.all$width))

###     tran all
PlotClusters(expr1=log10(d.all$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.all$s_ES_Som.Dicer_Rsc_tran_r2+1),
             ind = d.all$class,
             outpath=plot.path,
             name='O3.vs.S11.tran',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.all$width))


PlotClusters(expr1=log10(d.all$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.all$s_ES_WT_tran+1),
             ind = d.all$class,
             outpath=plot.path,
             name='O3.vs.WT.tran',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.all$width))


### utra 50
d.w = d.all[d.all$width >= 50,]


PlotClusters(expr1=log10(d.w$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.w$s_ES_Som.Dicer_Rsc_utra_r2+1),
             ind = d.w$class,
             outpath=plot.path,
             name='O3.vs.S11.utra',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.w$width))



PlotClusters(expr1=log10(d.w$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.w$s_ES_WT_utra+1),
             ind = d.w$class,
             outpath=plot.path,
             name='O3.vs.WT.utra',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.w$width))

###     tran w
PlotClusters(expr1=log10(d.w$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.w$s_ES_Som.Dicer_Rsc_tran_r2+1),
             ind = d.w$class,
             outpath=plot.path,
             name='O3.vs.S11.tran',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.w$width))


PlotClusters(expr1=log10(d.w$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.w$s_ES_WT_tran+1),
             ind = d.w$class,
             outpath=plot.path,
             name='O3.vs.WT.tran',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.w$width))






# ------------------------------------------------------------------- #
#without mirnas
d.mirna = d.all[d.all$class != 'miRNA',]
cols = c('red','black','darkgray')
#### utra all
PlotClusters(expr1=log10(d.mirna$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.mirna$s_ES_Som.Dicer_Rsc_utra_r2+1),
             ind = d.mirna$class,
             outpath=plot.path,
             name='O3.vs.S11.utra.womirna',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.mirna$width),
             cols = cols)



PlotClusters(expr1=log10(d.mirna$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.mirna$s_ES_WT_utra+1),
             ind = d.mirna$class,
             outpath=plot.path,
             name='O3.vs.WT.utra.womirna',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.mirna$width),
             cols = cols)

###     tran all
PlotClusters(expr1=log10(d.mirna$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.mirna$s_ES_Som.Dicer_Rsc_tran_r2+1),
             ind = d.mirna$class,
             outpath=plot.path,
             name='O3.vs.S11.tran.womirna',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.mirna$width),
             cols = cols)


PlotClusters(expr1=log10(d.mirna$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.mirna$s_ES_WT_tran+1),
             ind = d.mirna$class,
             outpath=plot.path,
             name='O3.vs.WT.tran.womirna',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.mirna$width),
             cols = cols)


### utra 50
d.w.mirna = d.mirna[d.mirna$width >= 50,]


PlotClusters(expr1=log10(d.w.mirna$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.w.mirna$s_ES_Som.Dicer_Rsc_utra_r2+1),
             ind = d.w.mirna$class,
             outpath=plot.path,
             name='O3.vs.S11.utra.womirna',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.w.mirna$width),
             cols = cols)



PlotClusters(expr1=log10(d.w.mirna$s_ES_MII.Dicer_Rsc_utra_r2+1),
             expr2=log10(d.w.mirna$s_ES_WT_utra+1),
             ind = d.w.mirna$class,
             outpath=plot.path,
             name='O3.vs.WT.utra.womirna',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.w.mirna$width),
             cols = cols)

###     tran w
PlotClusters(expr1=log10(d.w.mirna$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.w.mirna$s_ES_Som.Dicer_Rsc_tran_r2+1),
             ind = d.w.mirna$class,
             outpath=plot.path,
             name='O3.vs.S11.tran.womirna',
             xlab='DicerO-3',
             ylab='DicerS-11',
             mcw = min(d.w.mirna$width),
             cols = cols)


PlotClusters(expr1=log10(d.w.mirna$s_ES_MII.Dicer_Rsc_tran_r2+1),
             expr2=log10(d.w.mirna$s_ES_WT_tran+1),
             ind = d.w.mirna$class,
             outpath=plot.path,
             name='O3.vs.WT.tran.womirna',
             xlab='DicerO-3',
             ylab='DicerS-+/1',
             mcw = min(d.w.mirna$width),
             cols = cols)



# ------------------------------------------------------------------- #
# 3. Watanabe.data
library(rtracklayer)
lo.path='/common/USERS/vfranke/Base/Liftover/mm8/mm8ToMm9.over.chain'
chain = import.chain(lo.path)

wat.files = list.files(acc.path, recursive=T, full.names=T, patter='Wat')
wat.files = wat.files[str_detect(wat.files, 'endoSirna.txt')]
# wat = sapply(wat.files, function(x)BedToGRanges(read.table(x, header=T, sep='\t'), values=T))
names(wat) = str_replace(basename(wat.files),'.txt','')
watl = lapply(wat, function(x)liftOver(x, chain))
watgrl = lapply(watl, unlist)

watgr = unlist(GRangesList(watgrl))
names(watgr) = NULL
# d.wat = do.call(cbind, lapply(bams, function(x)GetOverlaps(watgr, x,'nh')))
d.wat = do.call(cbind, lapply(bams, function(x)countOverlaps(watgr, x,)))
total = unlist(lapply(l.mat, '[[', 'total'))
d.wat = t(t(d.wat)/total*1e6)
d.wat.out = data.frame(as.data.frame(watgr), d.wat)
write.table(d.wat.out, file.path(outpath, 'watanabe.rpm.txt'), row.names=F, col.names=T, quote=F, sep='\t')

wat.sel.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData/Watanabe_2008_Nature/WatanabeSelected.regions.txt'
wat.tab = read.table(wat.sel.path, header=T)
wat.sel = BedToGRanges(wat.tab)
s.wat = d.wat[countOverlaps(watgr, wat.sel) == 1,]
watgr.sel = watgr[countOverlaps(watgr, wat.sel) == 1]
ind1 = paste(as.character(seqnames(watgr.sel)), start(watgr.sel), end(watgr.sel))
ind2 = paste(wat.tab[,1], wat.tab[,2], wat.tab[,3])
m = match(ind2, ind1)
watout = data.frame(as.data.frame(watgr.sel),s.wat[,c(16,4,2,13,14,8:10,15,3,1,11,12,5:7)])
watout = watout[m,]
write.table(watout, file.path(outpath, 'wat.sel.txt'), row.names=F, col.names=T, quote=F, dec=',', sep='|')


# ------------------------------------------------------------------- #
# 3. Babiarz.data

library(rtracklayer)
lo.path='/common/USERS/vfranke/Base/Liftover/mm8/mm8ToMm9.over.chain'
chain = import.chain(lo.path)
chain = import.chain(lo.path)

bab.files = list.files(acc.path, recursive=T, full.names=T, patter='Bab')
bab.files = bab.files[str_detect(bab.files, 'txt')]
# bab = sapply(bab.files, function(x)BedToGRanges(read.table(x, header=T, sep='\t'), values=T))
names(bab) = str_replace(basename(bab.files),'.txt','')
babl = lapply(bab, function(x)liftOver(x, chain))
babgrl = lapply(babl, unlist)
for(i in names(babgrl)){
  print(i)
  values(babgrl[[i]])$class = i
}
babgr = unlist(GRangesList(babgrl))
names(babgr) = NULL
# d.bab = do.call(cbind, lapply(bams, function(x)GetOverlaps(babgr, x, 'nh')))
d.bab = do.call(cbind, lapply(bams, function(x)countOverlaps(babgr, x)))
total = unlist(lapply(l.mat, '[[', 'total'))
d.bab = t(t(d.bab)/total*1e6)
d.bab.out = data.frame(as.data.frame(babgr), d.bab)
write.table(d.bab.out, file.path(outpath, 'babanabe.rpm.txt'), row.names=F, col.names=T, quote=F, sep='\t')


### regions find out
regions = BedToGRanges(read.table(regions.path, header=T), values=T)
d.reg = do.call(cbind, lapply(bams, function(x)countOverlaps(regions, x)))
total = unlist(lapply(l.mat, '[[', 'total'))
d.reg.norm = data.frame(t(t(d.reg)/total*1e6))

fo = as.matrix(findOverlaps(regions, bams[['s_ES_MII.Dicer_Rsc_utra_r2']]))
fob = bams[['s_ES_MII.Dicer_Rsc_utra_r2']][fo[,2]]
fo.r = reduce(fob)
fo.br = as.matrix(findOverlaps(regions, fo.r))
fo.re = resize(fo.r, width=width(fo.r)+50)
foro = as.matrix(findOverlaps(fo.re, fo.re))
clu = ClusterGraphs(foro)
fod = do.call(cbind, lapply(bams, function(x)countOverlaps(fo.r, x)))
total = unlist(lapply(l.mat, '[[', 'total'))
fod.norm = data.frame(t(t(fod)/total*1e6))
aggregate(fod.norm, by=list(fo.br[,1]), FUN=sum)
aggregate(fod.norm, by=list(clu[,2]), FUN=sum)


svo.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData/SingleRegions/Svoboda.regions.txt'
svo.reg = BedToGRanges(read.table(svo.path, header=T))
cp = countOverlaps(svo.reg,bams[['s_ES_MII.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_MII.Dicer_Rsc_utra_r2']]) == '+'])
cm = countOverlaps(svo.reg,bams[['s_ES_MII.Dicer_Rsc_utra_r2']][strand(bams[['s_ES_MII.Dicer_Rsc_utra_r2']]) == '-'])
perc = cp/(cp+cm)
svod = data.frame(as.data.frame(svo.reg), plus=cp, minus=cm)
write.table(svod, file.path(outpath, 'custom.regions.strand.counts.txt'),row.names=F, col.names=T, quote=F, sep='\t')

s = getSeq(genome,svo.reg)
names(s) = paste(svod[,1], svod[,2], svod[,3], sep='_')
write.XStringSet(s, file.path(outpath, 'selected.reg.seq.fasta'))

reg = GRanges('chr17', IRanges(39980485, 39985676))
bams.single = lapply(bams, function(x)x[values(x)$nh == 1])
s = sapply(bams.single, function(x)sapply(c('+','-'), function(y)countOverlaps(reg, x[strand(x)==y])))
ts = t(s)
ts.norm = round(ts/total * 1e6,2)
write.table(ts, file.path(outpath, 'chr7.region.txt'), row.names=T, col.names=T, quote=F, sep=' ', dec=',')
write.table(ts, file.path(outpath, 'chr7.region.txt'), row.names=T, col.names=T, quote=F, sep=' ')
