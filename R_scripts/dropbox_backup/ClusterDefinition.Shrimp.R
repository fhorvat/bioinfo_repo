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
RpkmTab = function(tab, width){
  
  for(i in 1:ncol(tab)){
    tab[,i] = tab[,i]/(width*as.numeric(sum(tab[,i])))*10^6
  }
  return(tab)
}


#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES
inpaths='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/Mapped/ShrimpBest'

outpath= '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/ClusterDefinitionShrimp'
dir.create(outpath, showWarnings=F)

ens.genes.path = '/common/USERS/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/25012012.EnsemblGenesAll_BioMart.form.txt'

mirna.path = '/common/USERS/vfranke/Base/MirBase/mmu.19.gff3'

acc.path = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Data/AccessoryData'

#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
rpkm.cut = 10

genome.name = 'BSgenome.Mmusculus.UCSC.mm9'
#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES


# {{3}} MAIN
cuts = c(7,15,18, 20,23,34,35,49)
ncuts = c(0,1,4,10,20)

genome = GenomeLoader(genome.name)

bam.files = unlist(sapply(inpaths, list.files, pattern='bam$', full.names=T, recursive=T), use.names=F)
bam.files = bam.files[str_detect(bam.files, 'MII') | str_detect(bam.files, 'Som')]

seqlen = seqlengths(genome)
seqlen = seqlen[!str_detect(names(seqlen), '_')]
wins = MakeTillingWindows(100, seqlen)
regions = Reduce('c', wins)

l = list()
### saves the filtered bam files into a RData format
l.samps = list()
for(i in 1:length(bam.files)){
  
  bam.file = bam.files[i]
  name = BamName(bam.file)
  print(name)
  
  l.stats = list()
  # l.bam = list()
  # for(chr in chrs$chr){
  library(doMC)
  registerDoMC(2)
  l.bam = foreach(chr = names(seqlen))%dopar%{
    
    print(chr)
    l.counts = list()
    which = GRanges(chr, IRanges(1, seqlen[chr]))
    param = ScanBamParam(which=which, what=c('qname'))
    bam = readBamGappedAlignments(bam.file, param=param)
    
    d = data.table(qname=values(bam)$qname)
    n = d[,.N,by='qname']
    gr
  }
  r = Reduce('c', l.bam)
  
  l.stats$total = length(r)
  l.stats$total.unique = length(unique(names(r)))
  
  tab = table(names(r))
  NH = tab[match(names(r), names(tab))]
  w = unlist(width(r), use.names=F)
  wcuts = cut(w, c(min(w),cuts, max(w)), include.lowest=T)
  nhcuts = cut(NH, ncuts)
  cut.ind = paste(paste('w',wcuts,sep=''), paste('n',nhcuts,sep=''), sep=':')
  # sp = split(r, cut.ind)
  l.o = list()
  for(i in unique(cut.ind)){
    print(i)
    l.o[[i]] = countOverlaps(regions, r[cut.ind == i])
  }
  do = do.call(cbind, l.o)
  
  l.samps[[name]]$stats=l.stats
  l.samps[[name]]$counts = data.table(do)
}
# save(l, file=file.path(outpath, 'l.reads.RData'))
load(file.path(outpath, 'l.reads.RData'))
l.reads = lapply(l, '[[', 'reads')

# d.mirna = do.call(cbind, lapply(l.reads, function(x)countOverlaps(mirna, x)))
# df.mirna = cbind(as.data.frame(mirna), d.mirna)
# names(df.mirna)[1] = 'chr'
# write.table(df.mirna, file.path(outpath, 'dfmirna.txt'), row.names=F, col.names=T, quote=F, sep='\t')

# l.reads.new = endoapply(l.reads, function(x)x[countOverlaps(x, mirna) == 0])
# save(l, file=file.path(outpath, 'l.reads.new.RData'))

nam = lapply(l.reads.new, function(x)str_replace(names(x),'F3.+$','F3'))
tab = sapply(nam, table)
l.hits = lapply(names(tab), function(x)tab[[x]][match(nam[[x]], names(tab[[x]]))])
names(l.hits) = names(tab)


rangs = reduce(unlist(GRangesList(lapply(l.reads.new[-(1:4)], reduce))))
strand(rangs) = '*'
rangs = reduce(rangs)
l.reads.four = lapply(names(l.reads.new), function(x)l.reads.new[[x]][nam[[x]] %in% names(tab[[x]])[tab[[x]] <= 4]])
names(l.reads.four) = names(l.reads.new)
nam.four = lapply(l.reads.four, function(x)str_replace(names(x),'F3.+$','F3'))
d.all = do.call(cbind, lapply(names(l.reads.four), function(x){
  reads = l.reads.four[[x]]
  o = data.table(as.matrix(findOverlaps(rangs, reads)))
  name = nam.four[[x]][o$subjectHits]
  o$numHit = 1/l.hits[[x]][match(name, names(l.hits[[x]]))]
  h = o[, sum(numHit), by=queryHits]
  z = rep(0, length(rangs))
  z[h$queryHits] = h$V1
  z}))
d.all.rpkm = RpkmTab(d.all, width(rangs))
a = apply(d.all.rpkm, 1, max)
sel.ind = a > .5

rangs.sel = rangs[a > 1]
d.all.sel = d.all.rpkm[sel.ind,]
colnames(d.all.sel) = colnames(d.all)

# head(d.all[sel.ind,])[,1:5]
# head(d.one + d.two)[,1:5]
# head(d.one)[,1:5]
# head(d.two)[,1:5]

# 1 times
l.reads.one = lapply(names(l.reads.new), function(x)l.reads.new[[x]][nam[[x]] %in% names(tab[[x]][tab[[x]] == 1])])
names(l.reads.one) = names(l.reads.new)
d.one = do.call(cbind, lapply(l.reads.one, function(x)countOverlaps(rangs.sel, x)))
d.one.rpkm = RpkmTab(d.one, width(rangs.sel))

d.one.p = do.call(cbind, lapply(l.reads.one, function(x)countOverlaps(rangs.sel, x[strand(x) == '+'])))
d.one.m = do.call(cbind, lapply(l.reads.one, function(x)countOverlaps(rangs.sel, x[strand(x) == '-'])))

d.one.rat = log2((d.one.p+1) / (d.one.m +1))
d.one.rat = d.one.rat[,c(5:14)]

# 2 - 4 times
l.reads.two = lapply(names(l.reads.new), function(x)l.reads.new[[x]][nam[[x]] %in% names(tab[[x]][tab[[x]] > 1 & tab[[x]] <=4])])

names(l.reads.two) = names(l.reads.new)
nam.two = lapply(l.reads.two, function(x)str_replace(names(x),'F3.+$','F3'))
d.two = do.call(cbind, lapply(names(l.reads.two), function(x){
  reads = l.reads.two[[x]]
  o = data.table(as.matrix(findOverlaps(rangs.sel, reads)))
  name = nam.two[[x]][o$subjectHits]
  o$numHit = 1/l.hits[[x]][match(name, names(l.hits[[x]]))]
  h = o[, sum(numHit), by=queryHits]
  z = rep(0, length(rangs.sel))
  z[h$queryHits] = h$V1
  z}))

d.two.p = do.call(cbind, lapply(names(l.reads.two), function(x){
  reads = l.reads.two[[x]]
  reads = reads[strand(reads) == '+']
  o = data.table(as.matrix(findOverlaps(rangs.sel, reads)))
  name = nam.two[[x]][o$subjectHits]
  o$numHit = 1/l.hits[[x]][match(name, names(l.hits[[x]]))]
  h = o[, sum(numHit), by=queryHits]
  z = rep(0, length(rangs.sel))
  z[h$queryHits] = h$V1
  z}))

d.two.m = do.call(cbind, lapply(names(l.reads.two), function(x){
  reads = l.reads.two[[x]]
  reads = reads[strand(reads) == '-']
  o = data.table(as.matrix(findOverlaps(rangs.sel, reads)))
  name = nam.two[[x]][o$subjectHits]
  o$numHit = 1/l.hits[[x]][match(name, names(l.hits[[x]]))]
  h = o[, sum(numHit), by=queryHits]
  z = rep(0, length(rangs.sel))
  z[h$queryHits] = h$V1
  z}))

d.two.rat = log2((d.two.p+1)/(d.two.m+1))
d.two.rat = d.two.rat[,c(5:14)]

# library(DESeq)
# m = d.one[,c(5:14)]
# conds = str_replace(colnames(m),'_Rsc.+$','')
# cnt = newCountDataSet(m,conds)
# cnt = estimateSizeFactors(cnt)
# cnt = estimateDispersions(cnt)
# cnt = nbinomTest(cnt, "s_ES_MII.Dicer", "s_ES_Som.Dicer")

# more
l.reads.more = lapply(names(l.reads.new), function(x)l.reads.new[[x]][nam[[x]] %in% names(tab[[x]][tab[[x]] > 4])])
names(l.reads.more) = names(l.reads.new)

nam.more = lapply(l.reads.more, function(x)str_replace(names(x),'F3.+$','F3'))
d.more = do.call(cbind, lapply(names(l.reads.more), function(x){
  reads = l.reads.more[[x]]
  o = data.table(as.matrix(findOverlaps(rangs.sel, reads)))
  name = nam.more[[x]][o$subjectHits]
  o$numHit = 1/l.hits[[x]][match(name, names(l.hits[[x]]))]
  h = o[, sum(numHit), by=queryHits]
  z = rep(0, length(rangs.sel))
  z[h$queryHits] = h$V1
  z}))

r.more.filt = lapply(l.reads.more, function(x)x[countOverlaps(x, rangs) == 0])
r.more = reduce(unlist(GRangesList(lapply(r.more.filt[-(1:4)], reduce))))


# cluster.annotation
d = as.data.frame(rangs.sel, stringsAsFactors=F)
names(d)[1] = 'chr'
d$strand = as.character(d$strand)
d$WT.rsc.rat.one = rowMeans(d.one.rat[,4:7])
d$MII.rsc.rat.one = rowMeans(d.one.rat[,1:6])
d$WT.rsc.rat.two = rowMeans(d.two.rat[,4:7])
d$MII.rsc.rat.two = rowMeans(d.two.rat[,1:6])
pind = (rowSums(d[,(6:7)] >=.5) == 2) | (rowSums(d[,(8:9)] >=.5) == 2)
mind = (rowSums(d[,(6:7)] <= -.5) == 2) | (rowSums(d[,(8:9)] <= -.5) == 2)
d$strand[pind] = '+'
d$strand[mind] = '-'
d$strand[pind & mind] = 'b'


# cluster designation

l.g = list()
for(name in names(l.reads.two)){
  a = l.reads.two[[name]]
  n = nam.two[[name]]
  print(name)
  o = data.frame(as.matrix(findOverlaps(rangs.sel,a)))
  o$id = n[o[,2]]
  tb = table(o$id)
  o = o[o$id %in% names(tb)[tb>1],]
  o = o[,c(1,3)]
  s = split(o$id,o[,1])
  
  lo = list()
  for(i in 1:(length(s)-1)){
    n1 = names(s)[i]
    print(n1)
    for(j in (i+1):length(s)){
      n2 = names(s)[j]
      cat(n2,'\r')
      lo[[n1]][[n2]] = length(intersect(s[[n1]], s[[n2]]))
    }
  }
  l.g[[name]] = lo
}
save(l.g, file = file.path(outpath, 'l.g.RData'))
two.num = sapply(nam.two, function(x)length(unique(x)))
two.num = two.num/1e5
l.gc = lapply(names(l.g), function(x)lapply(l.g[[x]], function(y)y/two.num[x]))
names(l.gc) = names(l.g)
l.gc = lapply(l.gc, function(x)lapply(x, function(y)y[y>10]))
l.gg = lapply(l.gc, function(x){sind=sapply(x, length);x[sind!=0]})
a = l.gg[[1]]
l.dd = lapply(l.gg, function(x)lapply(x, function(y)data.frame(x2=names(y), count=y)))
l.dd = lapply(l.dd, function(x)data.frame(x1=rep(names(x), times=sapply(x,nrow)), do.call(rbind, x)))
l.dd = lapply(l.dd, apply, 2, as.numeric)
l.dd = lapply(l.dd, as.data.frame)

l.clust = list()
for(i in 1:length(l.sd)){
  l.s = l.dd[[i]][,1:2]
  v = vector()
  name = names(l.sd)[i]
  print(name)
  l.vec = list()
  inds = unique(l.s[,1])
  for(j in 1:length(inds)){
    
    ind = inds[j]
    cat(ind,'\r')
    if(ind %in% v){
      next
    }
    
    ng = l.s[l.s[,1] == ind,2]
    pg = ng
    while(! length(ng) == 0){
      
      for(k in ng){
        pg = c(pg, k)
        ng = unique(c(ng, unlist(l.s[l.s[,1] == k | l.s[,2] == k,])))
        ng = ng[ng != k]
        ng = ng[!ng %in% pg]
      }
      
    }
    v = unique(c(v, pg))
    l.vec[[j]]= unique(pg)
  }
  l.vec = l.vec[!sapply(l.vec, is.null)]
  l.clust[[name]] = l.vec
}
cd = lapply(l.clust, function(x)data.frame(cl=rep(1:length(x), times=sapply(x, length)),id=as.numeric(unlist(x))))
id = data.frame(id=1:nrow(d))
cd = lapply(names(cd), function(x){colnames(cd[[x]])[1] = x;cd[[x]]})
m = id
for(i in 1:length(cd)){
  print(i)
  m = merge(m, cd[[i]], by='id', all.x=T)
}

l.m = list()
for(i in 1:nrow(m)){
  l.c = list()
  cat(i,"\n")
  for(j in 2:ncol(m)){
    cat(j, '\r')
    l.c[[j]] = sum(m[,j] == m[i,j], na.rm=T)
  }
  l.m[[i]] = l.c
}
md = do.call(rbind, lapply(l.m, unlist))
tab = apply(md,1,table)
tab.num = apply(md,1, function(x)sum(x[-(1:4)] == 0))
d$connect = 'middle'
d$connect[tab.num >= 10] = 'None'
d$connect[tab.num <= 2] = 'Full'

tab.mean = round(apply(md, 1,function(x)mean(x[-(1:4)])),2)
d$conn.num = tab.mean

tab.two = apply(md,1, function(x)sum(x[-(1:4)] == 2))

l.t = list()
for(i in 2:ncol(m))
  l.t[[i]] = which(m[,i] == m[19,i])

m.all = list()
seen = vector()
for(i in 1:nrow(m)){
  cat(i, '\n')
  if(i %in% seen)
    next
  
  inds = vector()
  for(j in 7:ncol(m)){
    cat(i, '\r')
    inds = unique(c(inds, m$id[which(m[,j] == m[i,j])]))
  }
  m.all[[i]] = inds
  seen = c(seen, inds)
}
m.all = m.all[sapply(m.all, length) !=0]
d.all = data.frame(id = unlist(m.all), clust=rep(1:length(m.all), times=sapply(m.all, length)))
d$clust.num.all = 0
d$clust.num.all[d.all$id] = d.all$clust

getSeq(genome, GRanges('chr10', IRanges(20080833, 20080870)))
### annotation
annot = do.call(cbind, lapply(g.acc, function(x)ifelse(countOverlaps(rangs.sel, x) >0,1,0)))
annot.nam = apply(annot, 1, function(x)ifelse(sum(x==1) > 1, paste(colnames(annot)[x==1], collapse=':'), 'None'))
d$annot = annot.nam

d$class = 'NA'
d$class[d$strand %in% c('+', '-') & d$connect == 'None'] = 'miRNA'
d$class[!(d$strand %in% c('+', '-')) & d$connect == 'None'] = 'invRep'
d$class[(d$strand %in% c('+', '-')) & d$connect != 'None'] = 'endo.siRNA'


### mean expression

sets = c(1,1,2,2,3,3,3,3,3,3,4,4,4,4,5,5)
setnames = unique(str_replace(colnames(d.one), '_[tr|ut].*',''))
d.one.norm = t(t(d.one)/colSums(d.one)*1e5)
d.one.avg = t(apply(d.one.norm, 1, function(x)sapply(unique(sets), function(y)mean(x[sets==y]))))
colnames(d.one.avg) = paste(setnames,'.one.avg', sep='')
d.one.sd = t(apply(d.one.norm, 1, function(x)sapply(unique(sets), function(y)sd(x[sets==y]))))
colnames(d.one.sd) = paste(setnames,'.one.sd', sep='')

d.two.norm = d.two/two.num
d.two.avg = t(apply(d.two, 1, function(x)sapply(unique(sets), function(y)mean(x[sets==y]))))
colnames(d.two.avg) = paste(setnames,'.two.avg', sep='')
d.two.sd = t(apply(d.two, 1, function(x)sapply(unique(sets), function(y)sd(x[sets==y]))))
colnames(d.two.sd) = paste(setnames,'.two.sd', sep='')


### differential expression
library(DESeq)
samples = setnames[sets]
diff.ind=5:14
cnd.one = newCountDataSet(d.one[,diff.ind], samples[diff.ind])
cnd.one = estimateSizeFactors(cnd.one)
cnd.one = estimateDispersions(cnd.one, method='per-condition', sharingMode='fit-only')
test.one = nbinomTest(cnd.one, "s_ES_MII.Dicer_Rsc", "s_ES_Som.Dicer_Rsc")

test.one.d = test.one[,c(3,4,6,7)]
names(test.one.d) = c('mean.s_ES_MII.Dicer_Rsc.one','mean.s_ES_Som.Dicer_Rsc.two','log2FC.one','one.pval')
test.one.d$one.diff = 'None'
test.one.d$one.diff[which(test.one.d$one.pval <= .05 & test.one.d$log2FC.one > 0)] = 'Up'
test.one.d$one.diff[which(test.one.d$one.pval <= .05 & test.one.d$log2FC.one < 0)] = 'Down'
table(test.one.d$one.diff)


samples = setnames[sets]
d.two.count = apply(floor(d.two),2,as.integer)
cnd.two = newCountDataSet(d.two.count[,diff.ind], samples[diff.ind])
cnd.two = estimateSizeFactors(cnd.two)
cnd.two = estimateDispersions(cnd.two, method='per-condition', sharingMode='fit-only')
# fitType='locfit'
test.two = nbinomTest(cnd.two, "s_ES_MII.Dicer_Rsc", "s_ES_Som.Dicer_Rsc")

test.two.d = test.two[,c(3,4,6,7)]
names(test.two.d) = c('mean.s_ES_MII.Dicer_Rsc.two','mean.s_ES_Som.Dicer_Rsc.two','log2FC.two','two.pval')
test.two.d$two.diff = 'None'
test.two.d$two.diff[which(test.two.d$two.pval <= .05 & test.two.d$log2FC.two > 0)] = 'Up'
test.two.d$two.diff[which(test.two.d$two.pval <= .05 & test.two.d$log2FC.two < 0)] = 'Down'
table(test.two.d$two.diff)

dp = cbind(d, test.one.d, test.two.d)
cn=tapply(dp$chr, dp$clust.num.all, function(x)length(unique(x)))
dp$chr.num = cn[match(dp$clust.num.all, names(cn))]

plotDispEsts <- function( cds ){
  plot(
    rowMeans( counts( cds, normalized=TRUE ) ),
    fitInfo(cds)$perGeneDispEsts,
    pch = '.', log="xy" )
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}

library(Cairo)
CairoPNG(file.path(outpath, 'disp.est.one.only.png'), width=800, heigth=600)
plotDispEsts(cnd.one)
dev.off()
CairoPNG(file.path(outpath, 'disp.est.two.only.png'), width=800, heigth=600)
plotDispEsts(cnd.two)
dev.off()


### gene annotations
dg = BedToGRanges(cbind(dp[,1:3],id=1:nrow(dp)), values=T)
o = as.matrix(findOverlaps(dg, gr))
gene.ids = tapply(values(gr)$ens.gene.id[o[,2]], o[,1], function(x)paste(unique(x), collapse=':'))
gene.name = tapply(values(gr)$gen.name[o[,2]], o[,1], function(x)paste(unique(x), collapse=':'))
biotype = tapply(values(gr)$biotype[o[,2]], o[,1], function(x)paste(unique(x), collapse=':'))
dp$ens.gene.id = 'None'
dp$gene.name = 'None'
dp$biotype = 'None'
dp$ens.gene.id[as.numeric(names(gene.ids))] = gene.ids
dp$gene.name[as.numeric(names(gene.name))] = gene.name
dp$biotype[as.numeric(names(biotype))] = biotype

tss = gr

library(xlsx)

write.table(dp, file.path(outpath, 'cluster.annotation.txt'), col.names=T, row.names=T, quote=F, sep='\t')
write.table(d.one.avg, file.path(outpath, 'd.one.avg.txt'), col.names=T, row.names=T, quote=F, sep='\t')
write.table(d.one.sd, file.path(outpath, 'd.one.sd.txt'), col.names=T, row.names=T, quote=F, sep='\t')
write.table(d.two.avg, file.path(outpath, 'd.two.avg.txt'), col.names=T, row.names=T, quote=F, sep='\t')
write.table(d.two.sd, file.path(outpath, 'd.two.sd.txt'), col.names=T, row.names=T, quote=F, sep='\t')

xlsx.file = file.path(outpath, 'Dicer.xlsx')
options(java.parameters="-Xmx10024m")
library(xlsx)
gc()
write.xlsx(dp, file=xlsx.file, sheet='ClusterAnnotation', row.names=F, col.names=T)
write.xlsx(d.one.avg, file=xlsx.file, append=T, sheet='ExprOne.avg', row.names=F, col.names=T)
write.xlsx(d.one.sd, file=xlsx.file, append=T, sheet='ExprOne.sd', row.names=F, col.names=T)
write.xlsx(d.two.avg, file=xlsx.file, append=T, sheet='ExprTwo.avg', row.names=F, col.names=T)
write.xlsx(d.two.sd, file=xlsx.file, append=T, sheet='ExprTwo.sd', row.names=F, col.names=T)

table(dp$one.diff, dp$biotype)
rs = split(resize(dg, width=1, fix='center'), as.character(seqnames(dg)))
dists= lapply(rs, function(x){
  distm=as.matrix(dist(start(x)))
  distm[col(distm) == row(distm)] = Inf
  dists = apply(distm, 1, min)
  dists})
dists = unlist(dists)
CairoPNG(file.path(outpath, 'dists.png'), width=800, height=600)
hist(log10(dists), breaks=50, col='darkorange')
dev.off()

### annotation counts
annot.regs = reduce(unlist(GRangesList(lapply(g.acc, function(x){values(x)=NULL;x}))))
annot.d = do.call(cbind, lapply(g.acc, function(x)ifelse(countOverlaps(annot.regs, x) >0,1,0)))
annot.d = apply(annot, 1, function(x)ifelse(sum(x==1) > 1, paste(colnames(annot)[x==1], collapse=':'), 'None'))
annot.d.one  = do.call(cbind, lapply(l.reads.one, function(x)countOverlaps(annot.regs, x)))
r = d.one.sd/d.one.avg


annot.d.two = do.call(cbind, lapply(names(l.reads.two), function(x){
  reads = l.reads.two[[x]]
  o = data.table(as.matrix(findOverlaps(annot.regs, reads)))
  name = nam.two[[x]][o$subjectHits]
  o$numHit = 1/l.hits[[x]][match(name, names(l.hits[[x]]))]
  h = o[, sum(numHit), by=queryHits]
  z = rep(0, length(annot.regs))
  z[h$queryHits] = h$V1
  z}))
annot.ind = rowSums(annot.d.one) != 0 | rowSums(annot.d.two) != 0
annot.d.one = annot.d.one[annot.ind,]
annot.d.two = annot.d.two[annot.ind,]
annot.regs = annot.regs[annot.ind]

#/{{3}} MAIN
#/{2} CODE

save.image(file=file.path(outpath,'workspace.RData'))
load(file=file.path(outpath,'workspace.RData'))

path = file.path(outpath, 'cluster.annotation.txt')
d = read.table(path, header=T)
d.ran = BedToGRanges(cbind(d[,1:3], id = 1:nrow(d)), values=T)
d.ran = d.ran[order(values(d.ran)$id)]

tss = gr
sind = as.character(strand(tss)) == '+'
tss[sind] = resize(tss[sind], fix = 'start', width=1000)
tss[!sind] = resize(tss[!sind], fix = 'end', width=1000)
tss = ifelse(countOverlaps(d.ran, tss) >0, 'Prom', 'None')
body = ifelse(countOverlaps(d.ran, gr) >0, 'Body', 'None')
d$position = 'None'
d$position[body == 'Body'] = 'Body'
d$position[tss == 'Prom'] = 'Prom'

gs = gr[values(gr)$ens.gene.id  %in% d$ens.gene.id]
gs = as.data.frame(gs)
d$gene.strand = gs$strand[match(d$ens.gene.id, gs$ens.gene.id)]
d$gene.strand[is.na(d$gene.strand)] = '*'
write.table(d, file.path(outpath, 'cluster.annotation.new.txt'), row.names=F, col.names=T, quote=F, sep='\t', dec=',')
