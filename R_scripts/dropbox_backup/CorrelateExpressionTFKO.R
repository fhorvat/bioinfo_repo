### INFO: R Script
### DATE: 30.05.2014
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
library(affy)
library(gcrma)

#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS
#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES

# {{{1}}} PATH VARIABLES

TF.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/AccessoryData/MicroarrayData'

mt.id.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/MT_Transcriptome/MT.EnsGeneID.moe430.2.txt'
mt.id.best.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/ReceivedData/Petr/MTTranscriptionFactor/MTSelGenesOOcyte.Dev.140530.txt'

mt.genes.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/MT_Transcriptome/MT.one.star.txt'

outpath = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/MT_Transcriptome/MT.TFknockouts'

#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS

#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES


# {{3}} MAIN

TF.files = list.files(TF.path, full.names=TRUE, recursive=TRUE, pattern='KO')
mt = read.table(mt.id.path, header=TRUE, sep='\t')
mt = unique(mt)

mt.id = read.table(mt.id.best.path, header=F, sep='\t')

# -------------------------------------------- #
# LHX8 - Nobox
file = TF.files[str_detect(TF.files, 'Nobox')]
workdir = dirname(file)
Assigner(file, 'g')

pdata = do.call(rbind, lapply(g, pData))
links = pdata$supplementary_file
registerDoMC(length(links))
foreach(i = 1:length(links))%dopar%{
  
  link = links[i]
  outfile = file.path(dirname(file), basename(link))
  download.file(link, outfile)
}
cel.files = list.celfiles(workdir, pattern='GSM')
pdata = pdata[basename(links) %in% basename(cel.files),]
rownames(pdata) = basename(cel.files)
pdata$title = c('Nobox.KO','Nobox.KO','LHX8.KO','LHX8.KO','LHX8.KO','WT','WT')
affy = ReadAffy(filenames=cel.files, phenoData=pdata)

norm = gcrma(affy)
np = pData(norm)

nob = exprs(norm)
colnames(nob) = np$title[match(colnames(nob), rownames(np))]

nob.ind = str_detect(colnames(nob), 'Nobox')
wt.ind = str_detect(colnames(nob), 'WT')
M.nob = rowMeans(nob[,nob.ind]) - rowMeans(nob[,wt.ind])
A.nob = (rowMeans(nob[,nob.ind]) + rowMeans(nob[,wt.ind]))/2

lhx.ind = str_detect(colnames(nob), 'LHX')
M.lhx = rowMeans(nob[,lhx.ind]) - rowMeans(nob[,wt.ind])
A.lhx = (rowMeans(nob[,lhx.ind]) + rowMeans(nob[,wt.ind]))/2

sind = rownames(nob) %in% mt[,2] | rownames(nob) %in% mt[,3]
bind = rownames(nob) %in% mt.id[mt.id[,1]=="MT-driven",3]
pdf(file.path(outpath, 'Nobox.LHX8.KO.pdf'), width=10, height=8)
plot(A.nob,M.nob, pch=20, cex=.5, xlab='log2( NoboxKO + WT ) / 2', ylab='log2( NoboxKO / WT )', col='gray', main='NoboxKO')
points(A.nob[sind],M.nob[sind], pch=20, cex=1, col='red')
points(A.nob[bind],M.nob[bind], pch=20, cex=1, col='blue')
abline(h=0, lwd=2)

plot(A.lhx, M.lhx, pch=20, cex=.5, xlab='log2( LHX8.KO + WT ) / 2', ylab='log2( LHX8.KO / WT )',col='gray', main='LHX8 KO')
points(A.lhx[sind],M.lhx[sind], pch=20, cex=.75, col='red')
points(A.lhx[bind],M.lhx[bind], pch=20, cex=1, col='blue')
abline(h=0, lwd=2)

plot(M.nob, M.lhx, pch=20, cex=.5, xlab='Nobox KO / WT', ylab=' LHX8 KO / WT ', col='gray', main='Nobox VS LHX8')
points(M.nob[sind],M.lhx[sind], pch=20, cex=1, col='red')
points(M.nob[bind],M.lhx[bind], pch=20, cex=1, col='blue')
abline(h=0, lwd=1)
abline(v=0, lwd=1)
dev.off()





# -------------------------------------------- #
# Shlh
mt.ill.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/MT_Transcriptome/MT.EnsGeneID.ILL.MouseWG6v2.txt'
mt.ill = read.table(mt.ill.path, header=T, sep='\t')
mt.id = read.table(mt.id.path, sep='\t', header=TRUE)
mt.id = mt.id[mt.id$type == 'MT-driven',]
mt.ill.id = mt.ill[mt.ill[,1] %in% mt.id$ens.gene.id,]

file = TF.files[str_detect(TF.files, 'Shlh')]
workdir = dirname(file)
Assigner(file, 'g')

pdata = do.call(rbind, lapply(g, pData))
shl = exprs(g[[1]])
title = str_replace(pdata$title, '-.+', '')
title = str_replace(title, '\\s+', '.')
colnames(shl) = title[match(colnames(shl), str_replace(rownames(pdata),'^.+\\.',''))]

sind1  = str_detect(colnames(shl), 'Sohlh1')
wind = str_detect(colnames(shl), 'Wildtype')
M.shl1 = rowMeans(shl[,sind1]) - rowMeans(shl[,wind])
A.shl1 = (rowMeans(shl[,sind1]) + rowMeans(shl[,wind]))/2

sind2 = str_detect(colnames(shl), 'Sohlh2')
M.shl2 =  rowMeans(shl[,sind2]) - rowMeans(shl[,wind])
A.shl2 = (rowMeans(shl[,sind2]) + rowMeans(shl[,wind]))/2

sind = rownames(shl) %in% mt.ill[,2]
bind = rownames(shl) %in% mt.ill.id[,2]
pdf(file.path(outpath, 'Sohlh.KO.pdf'), width=10, height=8)
plot(A.shl1, M.shl1, pch=20, cex=.5, xlab='log2( Sohlh1.KO + WT ) / 2', ylab='log2( Sohlh1.KO / WT )', col='gray', main='Sohlh1.KO')
points(A.shl1[sind],M.shl1[sind], pch=20, cex=1, col='red')
points(A.shl1[bind],M.shl1[bind], pch=20, cex=1, col='blue')
abline(h=0, lwd=2)

plot(A.shl2, M.shl2, pch=20, cex=.5, xlab='log2( Sohlh2.KO + WT ) / 2', ylab='log2( Sohlh2.KO / WT )', col='gray', main='Sohlh2.KO')
points(A.shl2[sind],M.shl2[sind], pch=20, cex=1, col='red')
points(A.shl2[bind],M.shl2[bind], pch=20, cex=1, col='blue')
abline(h=0, lwd=2)
dev.off()

# -------------------------------------------- #
# FoxO3

mt.ill.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/MT_Transcriptome/MT.EnsGeneID.ILL.MouseWG6v2.txt'
mt.ill = read.table(mt.ill.path, header=T, sep='\t')
mt.id = read.table(mt.id.path, sep='\t', header=TRUE)
mt.id = mt.id[mt.id$type == 'MT-driven',]
mt.ill.id = mt.ill[mt.ill[,1] %in% mt.id$ens.gene.id,]

file = TF.files[str_detect(TF.files, 'FoxO3')]
workdir = dirname(file)
Assigner(file, 'g')
pdata = do.call(rbind, lapply(g, pData))
fox = exprs(g[[1]])
fox = log2(fox)

title = str_replace(pdata$title, 'Ovary-', '')
colnames(fox) = title[match(colnames(fox), str_replace(rownames(pdata),'^.+\\.',''))]
find = str_detect(colnames(fox), 'NoTg-Foxo3-/-')
wind = str_detect(colnames(fox), 'NoTg-Foxo3\\+/-')
M.fox = rowMeans(fox[,find]) - rowMeans(fox[,wind])
A.fox = (rowMeans(fox[,find]) + rowMeans(fox[,wind]))/2

sind = rownames(fox) %in% mt.ill[,2]
bind = rownames(fox) %in% mt.ill.id[,2]
pdf(file.path(outpath, 'Foxo3.KO.pdf'), width=10, height=8)
plot(A.fox, M.fox, pch=20, cex=.5, xlab='log2( Foxo3.KO + WT ) / 2', ylab='log2( Foxo3.KO / WT )', col='gray', main='Foxo3.KO')

points(A.fox[sind],M.fox[sind], pch=20, cex=1, col='red')
points(A.fox[bind],M.fox[bind], pch=20, cex=1, col='blue')
abline(h=0, lwd=2)
dev.off()

# -------------------------------------------- #
# Figla


# -------------------------------------------- #
# Shlh double ko

mt.ill.path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Results/MT_Transcriptome/MT.EnsGeneID.ILL.MouseWG6v2.txt'
mt.ill = read.table(mt.ill.path, header=T, sep='\t')
mt.id = read.table(mt.id.path, sep='\t', header=TRUE)
mt.id = mt.id[mt.id$type == 'MT-driven',]
mt.ill.id = mt.ill[mt.ill[,1] %in% mt.id$ens.gene.id,]

path = '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/AccessoryData/MicroarrayData/GSE24815_Shlh_DoubleKO'
id = str_replace(basename(path),'_.+','')
print(id)
library(GEOquery)
g = getGEO(id)
g = g[[1]]
p = pData(g)
snames = str_replace(p$title, '-.+', '')
shl  = exprs(g)
shl = log2(shl)
colnames(shl) = snames[match(p$geo_accession, colnames(shl))]

sind1  = str_detect(colnames(shl), 'DKO')
wind = str_detect(colnames(shl), 'Wildtype')
M.shl1 = rowMeans(shl[,sind1]) - rowMeans(shl[,wind])
A.shl1 = (rowMeans(shl[,sind1]) + rowMeans(shl[,wind]))/2


sind = rownames(shl) %in% mt.ill[,2]
bind = rownames(shl) %in% mt.ill.id[,2]
pdf(file.path(outpath, 'Sohlh.DoubleKO.pdf'), width=10, height=8)
plot(A.shl1, M.shl1, pch=20, cex=.5, xlab='log2( Sohlh1.Double.KO + WT ) / 2', ylab='log2( Sohlh1.KO / WT )', col='gray', main='Sohlh1.KO')
points(A.shl1[sind],M.shl1[sind], pch=20, cex=1, col='red')
points(A.shl1[bind],M.shl1[bind], pch=20, cex=1, col='blue')
abline(h=0, lwd=2)
dev.off()

#/{{3}} MAIN
#/{2} CODE
