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
source(file.path(lib.path, 'ScanLib.R'))

library(Rsamtools)
library(doMC)
library(data.table)
#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS
# ------------------------------------------------------------------------ #
WidthDist = function(cwall, x, outpath, name, axis.lwd=3, lwd=6, cols, main='All', legend, ylim.shift=0){
  
  x = as.vector(x)
  mi = min(x)
  mx = max(x)
  outfile = file.path(outpath, paste(name,'.lw',lwd, '.alw',axis.lwd,'.mi',mi,'.mx',mx,',.pdf', sep=''))
  CairoPDF(outfile, width=9, height=6)
  xlim = range(x)
  ylim=c(0, max(cwall+ylim.shift))
  for(i in 1:ncol(cwall)){
    
    print(i)
    v = as.vector(cwall[,i])
    if(i == 1){
      plot(x, v,  type='l', ylim = ylim, xlim=xlim, col=cols[i], xlab='Read width', ylab='Frequency', main=main, lwd=lwd, axes=F)
      axis(1, xlim[1]:xlim[2], cex.axis=1, lwd=axis.lwd)
      axis(2,seq(ylim[1],ylim[2],0.02), cex.axis=1, lwd=axis.lwd)
    }else{
      lines(x, v, col=cols[i], lwd=lwd)
    }
    
  }
  legend('topright', legend = legend, fill=cols, bty='n')
  dev.off()
}


# ------------------------------------------------------------------------ #
GetStats = function(x){
  
  cat('Doing the stats...\n')
  stats = list()
  stats$total = length(x$qname)
  stats$total.uniq = x[,length(unique(qname))]
  stats$width.dist.all = x[,.N, by=width]
  stats$annot.dist.all = x[,.N, by=ind]
  stats$annot.dist.frac = x[,round(sum(1/nh),2), by=ind]
  stats$mult.all = x[, .N, by=nh]
  stats$mult = stats$mult.all[,uniq := N/nh]
  
  dx = x[,unique(paste(qname, ind))]
  dx = unlist(strsplit(dx, ' '))
  dx.tab = table(dx[1:length(dx) %% 2 == 0])
  stats$annot.dist.uniq = data.table(class=names(dx.tab), counts=dx.tab)
  
  wx = x[,unique(paste(qname, width))]
  wx = unlist(strsplit(wx, ' '))
  wx.tab = table(wx[1:length(wx) %% 2 == 0])
  stats$width.uniq = data.table(width=names(wx.tab), counts=wx.tab)
  cat('Returning the list...\n')
  return(stats)
}

#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES
inpaths='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/Mapped/ShrimpStrata/'

cut.path='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/Cutadapt_Trimmed2'

outpath='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/MappingStats/ShrimpStrata/2012-10-03_Stats'
# outpath = file.path(outdir, OutpathDate('Stats'))
dir.create(outpath, showWarnings=F)

rdata=file.path(outpath, 'RData4')
dir.create(rdata, showWarnings=F)

wpath=file.path(outpath, 'Width')
dir.create(wpath, showWarnings=F)

apath=file.path(outpath, 'Annotation')
dir.create(apath, showWarnings=F)

stats.rdata.path = file.path(outpath, 'RDataStats')
dir.create(stats.rdata.path, showWarnings=F)


ens.genes.path = '/common/USERS/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/25012012.EnsemblGenesAll_BioMart.form.txt'
ens.genes.gtf.path = '/common/USERS/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/22012012_EnsemblGenes_BioMart.form.gtf'

mirna.path = '/common/USERS/vfranke/Base/MirBase/mm9.LO.mm10.mirbase.form.gff'

acc.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData'

repeat.path='/common/USERS/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.repeatmasker.RData'

sampnames.path='/home/members/vfranke/Projects/Dicer_RNASeq_29062012/Scripts/SampleLabels.R'

cutadapt.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/Cutadapt_Trimmed'

exons.path = '/common/USERS/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/22012012_EnsemblGenes_BioMart.form.gtf'

#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
registerDoMC(6)
# registerDoMC(9)
#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES


# {{3}} MAIN
source(sampnames.path)

cuts = c(7,15,18, 20,23,34,35,49)

bam.files = list.files(inpaths, pattern='bam$', full.names=T, recursive=T)

mirna = read.table(mirna.path, header=F)
mirna[,9] = str_replace_all(mirna[,9],'=',' ')
mirna = mirna[!(is.na(mirna[,4]) | is.na(mirna[,5])),]

mirna = GffToGRanges(mirna)

exons = read.table(exons.path, header=F, sep='\t')
exons = GffToGRanges(exons)

genes = read.table(ens.genes.path, header=T)
genes = BedToGRanges(genes, values=T)
genes = genes[!values(genes)$biotype %in% c('pseudogene','miRNA')]
exons = exons[values(exons)$gene_id %in% values(genes)$ens.gene.id]

gr = reduce(genes)
ex = reduce(exons)

Assigner(repeat.path, 'reps')
reps = updateObject(reps)
reps = reps[!str_detect(as.vector(seqnames(reps)), 'random')]
seqlevels(reps) = unique(as.character(seqnames(reps)))
reps = reduce(reps)

l.annot = list(
  miRNA=mirna,
  Genes=ex,
  Repeats=reps)

rm(reps, gr, genes, exons, ex);gc()

for(i in 1:length(bam.files)){
  
  bam.file = bam.files[i]
  name=BamName(bam.file)
  print(name)
  chrs = chrFinder(bam.file)
  # chrs = chrs[str_detect(chrs$chr, 'chr'),]
  
  if(str_detect(bam.file, 'ES_Som')){
    registerDoMC(4)
  }else{
    registerDoMC(6)
  }
  
  l = foreach(chr = chrs$chr)%dopar%{
    
    print(chr)
    cat('Reading in the files...\n')
    which.ranges = GRanges(chr, IRanges(1, chrs$chr.len[chrs$chr == chr]))
    param = ScanBamParam(what=c('qname'),which=which.ranges, flag=scanBamFlag(isUnmappedQuery=FALSE), tag=c('NM'))
    bam = readGappedAlignments(bam.file, param = param)
    bam = bam[values(bam)$NM <= 2]
    
    cat('Overlapping annotations...\n')
    d = do.call(cbind, lapply(l.annot, function(x)suppressWarnings(countOverlaps(bam, x))))
    d[d>1] = 1
    d = t(t(d) * 1:ncol(d))
    
    cat('Annotating Reads...\n')
    d[d == 0] = max(d)+1
    da = t(data.frame(d))
    mi = do.call(pmin, split(da, 1:ncol(d)))
    ind = c(names(l.annot), 'Other')[mi]
    
    dt = data.table(width = qwidth(bam), qname=values(bam)$qname, ind=ind, NM=values(bam)$NM)
    rm(ws,s, bam, d, da, s, mi, ind);gc()
    
    return(dt)
  }
  cat('Rbindlist...\n')
  r = rbindlist(l)
  
  cat('Summarizing the data...\n')
  
  cat('Multimapping...\n')
  tab = r[,.N, by=qname]
  m = match(r$qname, tab$qname)
  r$nh = tab$N[m]
  r$cut = cut(r$width, cuts, include.lowest=T)
  
  
  cat('Saving\n')
  save(r, file=file.path(rdata, paste(name, 'data.table.RData', sep='.')))
  
  rm(r, tab, m);gc()
}
# number of mapped reads


# ------------------------------------------------------------------------ #
# Statistics RData
rdata.files = list.files(rdata, full.names=T, pattern="RData")
rdata.files = rdata.files[str_detect(rdata.files, 's_')]

registerDoMC(9)
foreach(i = 1:length(rdata.files))%dopar%{
  
  rdata.file = rdata.files[i]
  name = str_replace(basename(rdata.file),'.filt.data.table.RData','')
  print(name)
  Assigner(rdata.file, 'l')
  
  jin =  sort(l[,unique(NM)])[1:3]
  stats.all = list()
  for(j in jin){
    print(j)
    jind = as.character(j)
    
    cat('Taking all reads...\n')
    lc = subset(l, subset = NM <= j)
    tablc = lc[,.N, by=qname]
    m = match(lc$qname, tablc$qname)
    lc$nh = tablc$N[m]
    stats.all[[paste(jind, 'all', 'all',sep='.')]] = GetStats(lc)
    
    cat('Taking the single mapping reads...\n')
    lc.all.1 = subset(lc, subset = nh ==1)
    stats.all[[paste(jind, 'all', '1',sep='.')]] = GetStats(lc.all.1)
    
    cat('Taking the 21 - 23 reads...\n')
    lc21 = subset(lc, subset = cut == '(20,23]')
    tablc21 = lc21[,.N, by=qname]
    m = match(lc21$qname, tablc21$qname)
    lc21$nh = tablc21$N[m]
    stats.all[[paste(jind, '21','all', sep='.')]] = GetStats(lc21)
    
    cat('Taking the 21 - 23 readsm multimapping 2...\n')
    lc21.2 = subset(lc21, subset = nh <= 2)
    stats.all[[paste(jind, '21','2', sep='.')]] = GetStats(lc21.2)
    
    cat('Taking the 21 - 23 readsm multimapping 1...\n')
    lc21.1 = subset(lc21, subset = nh == 1)
    stats.all[[paste(jind, '21','1', sep='.')]] = GetStats(lc21.1)
    
    rm(list=ls()[grepl('lc',ls())]);gc()
  }
  
  save(stats.all, file = file.path(stats.rdata.path, paste(name,'stats','RData', sep='.')))
}



# ------------------------------------------------------------------------ #
stats.files = list.files(stats.rdata.path, full.names=T, pattern='.RData')
stats.files = stats.files[str_detect(stats.files, 's_')]
l.d = list()
for(i in 1:length(stats.files)){
  
  stat.file = stats.files[i]
  name = str_replace(basename(stat.file),'.stats.RData','')
  name = str_replace(basename(name),'.data.table.RData','')
  print(name)
  Assigner(stat.file, 'l')
  
  a = as.data.frame(Reduce(function(x,y)merge(x,y,by='ind',all=T), lapply(l, '[[', 'annot.dist.frac')))
  colnames(a)[-1] = names(l)
  
  w = as.data.frame(Reduce(function(x,y)merge(x,y,by='width', all=T), lapply(l, '[[', 'width.uniq')))
  w[is.na(w)] = 0
  colnames(w)[-1] = names(l)
  w = w[,-c(4:6,9:11,14:16)]
  
  l.d[[name]]$counts = a
  l.d[[name]]$widths = w
}


### total counts
total = unlist(lapply(l.stats.sub, '[[', 'total'))
total.uniq = unlist(lapply(l.stats.sub, '[[', 'total.uniq'))
total.one = sapply(l.stats.sub, function(x){a=data.frame(x[['mult.all']]);a[a$N==1,2]})

# ----------------------------------------------------------------------------------------- #
# WIDTH DISTRIBUTION


# ------------------------------------------------------ #
### width distribution all
widths = lapply(l.d, '[[', 'widths')
widths = lapply(widths, function(x)as.matrix(x[4:15,]))

ncols = ncol(widths[[1]])
wutra = widths[str_detect(names(widths), 'utra')]
wutra = wutra[-c(3,5,6)]
wutra = wutra[c(5,2,1,3,4)]
width = as.numeric(wutra[[1]][,1])
for(i in 2:ncols){
  
  setname = colnames(wutra[[1]])[i]
  wd = do.call(cbind, lapply(wutra, function(x)as.numeric(x[,i])))
  wd = t(t(wd)/colSums(wd))
  legend = sapply(sampname[colnames(wd)], '[[', 'expr')
  cols = c('cornflowerblue', 'darkorange', 'red', 'black', 'gray')
  WidthDist(wd, width, cols=cols, main='Utra', outpath=wpath, name=paste('Utra.Thick', setname, sep='.'), legend=legend, ylim.shift=.02)
  WidthDist(wd, width, cols=cols, main='Utra', outpath=wpath, name=paste('Utra.Thin', setname, sep='.'), legend=legend, lw=3, ylim.shift=.02, axis.lwd=1.5)
  
}


ncols = ncol(widths[[1]])
wtran = widths[str_detect(names(widths), 'tran')]
wtran = wtran[-c(3,5,6)]
wtran = wtran[c(5,2,1,3,4)]
width = as.numeric(wtran[[1]][,1])
for(i in 2:ncols){
  
  setname = colnames(wtran[[1]])[i]
  wd = do.call(cbind, lapply(wtran, function(x)as.numeric(x[,i])))
  wd = t(t(wd)/colSums(wd))
  legend = sapply(sampname[colnames(wd)], '[[', 'expr')
  cols = c('cornflowerblue', 'darkorange', 'red', 'black', 'gray')
  WidthDist(wd, width, cols=cols, main='Tran', outpath=wpath, name=paste('Tran.Thick', setname, sep='.'), legend=legend, ylim.shift=.02)
  WidthDist(wd, width, cols=cols, main='Tran', outpath=wpath, name=paste('Tran.Thin', setname, sep='.'), legend=legend, lw=3, ylim.shift=.02, axis.lwd=1.5)
  
}


# ------------------------------------------------------ #


# ----------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------- #           # ANNOTATION DISTRIBUTION


# ------------------------------------------------------ #
# Mult Annotation
# annot = lapply(l.stats.sub, '[[', "annot.dist.all")
annot = lapply(l.d, function(x)x$counts[,c(1,which(str_detect(colnames(x$counts),'21')))])

ncols = ncol(annot[[1]])
annot.tran = annot[str_detect(names(annot), 'tran')]
annot.tran = annot.tran[c(8,2,1,6,7,3,4,5)]
annot.name = l.d[[1]]$counts[,1]
sort.ind = c(2,1,4,3)
for(i in 2:ncols){
  
  setname = colnames(annot.tran[[1]])[i]
  print(setname)
  dat = do.call(cbind, lapply(annot.tran, '[[', i))
  rownames(dat) = annot.name
  dat = dat[sort.ind,]
  total = colSums(dat)
  dat = t(t(dat)/total)
  
  write.table(dat, file.path(apath, paste(setname,'.tran.txt',sep='.')), row.names=T, col.names=T, quote=F, dec=',', sep=' ')
  names.arg = sapply(sampname[colnames(dat)], '[[', 'expr')
  annot.cols = c('black', gray(c(.3,.6)),'white')
  bar.lwd=3.5
  CairoPNG(file.path(apath, paste('Dicer.21-23.Trans',setname,'png', sep='.')),width=900, height=800)
  layout(matrix(1:2, nrow=2), heights = c(4,1))
  par(lwd=bar.lwd)
  barplot(dat, names.arg=names.arg, ylab='Frequency', main='Transfected Annotation', , col=annot.cols, mar=rep(0,4), oma=rep(0,4),lwd=bar.lwd)
  axis(1, pos=0.00005, labels=F, at=c(-.15,.71,c(0.71+cumsum(rep(1.2,7)))), lwd=bar.lwd)
  plot(1,1,axes=F, xlab='', ylab='', col='white')
  par(lwd=1)
  legend('center', legend = as.character(rownames(dat)), bty='n', fill=annot.cols, cex=1, horiz=T, text.width=.1, pt.cex=1.5, inset=c(10,0))
  dev.off()
}
# annot = t(t(annot)/total)
# means tran


ncols = ncol(annot[[1]])
annot.utra = annot[str_detect(names(annot), 'utra')]
annot.utra = annot.utra[c(8,2,1,6,7,3,4,5)]
annot.name = l.d[[1]]$counts[,1]
sort.ind = c(2,1,4,3)
for(i in 2:ncols){
  
  setname = colnames(annot.utra[[1]])[i]
  print(setname)
  dat = do.call(cbind, lapply(annot.utra, '[[', i))
  rownames(dat) = annot.name
  dat = dat[sort.ind,]
  total = colSums(dat)
  dat = t(t(dat)/total)
  
  write.table(dat, file.path(apath, paste(setname,'utra.txt',sep='.')), row.names=T, col.names=T, quote=F, dec=',', sep=' ')
  names.arg =sapply(sampname[colnames(dat)], '[[', 'expr')
  annot.cols = c('black', gray(c(.3,.6)),'white')
  bar.lwd=3.5
  CairoPNG(file.path(apath, paste('Dicer.21-23.Utra',setname,'png',sep='.')),width=900, height=800)
  layout(matrix(1:2, nrow=2), heights = c(4,1))
  par(lwd=bar.lwd)
  barplot(dat, names.arg=names.arg, ylab='Frequency', main='Untransfected Annotation', , col=annot.cols, mar=rep(0,4), oma=rep(0,4),lwd=bar.lwd)
  axis(1, pos=0.00005, labels=F, at=c(-.15,.71,c(0.71+cumsum(rep(1.2,7)))), lwd=bar.lwd)
  plot(1,1,axes=F, xlab='', ylab='', col='white')
  par(lwd=1)
  legend('center', legend = as.character(rownames(dat)), bty='n', fill=annot.cols, cex=1, horiz=T, text.width=.1, pt.cex=1.5, inset=c(10,0))
  dev.off()
}


# ------------------------------------------------------ #




# ----------------------------------------------------------------------------------------- #










#/{{3}} MAIN
#/{2} CODE
