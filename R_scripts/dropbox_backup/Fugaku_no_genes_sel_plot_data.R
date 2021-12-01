#!/common/WORK/vfranke/bin/R/R-3.1.0/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/")

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
library(ggplot2)
library(readr)

#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS

#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES 

# {{{1}}} PATH VARIABLES
inpath='/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/R_objects/no_genes_sel'

outpath='/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/results'

# samp.ord = '/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/vfranke_documentation/SampleOrder.R'

tot.path = '/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/vfranke_documentation/TotCounts.R'


#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
registerDoMC(21)

source(tot.path)

source(samp.ord)
#/{{{2}}} SCRIPT PARAMS

#/{{2}} INPUT VARIABLES


# {{3}} MAIN


rfiles = list.files(inpath, full.names=TRUE, pattern='RData')
l = foreach(i = 1:length(rfiles))%dopar%{
  
  print(i)
  Assigner(rfiles[i], 'a')
  return(a)
}
names(l) = str_replace(basename(rfiles),'.SpliceDonAcc.RData','')
# l = l[ord]

sl = lapply(l, function(x)x[[2]])
sl = lapply(names(sl), function(x){a=sl[[x]];a$samp=x;a})
sl = rbindlist(sl)
sl$samp = factor(sl$samp, levels=c("s_GV.WE", "s_MII.WE", "s_1cell.WE", "s_1cell.WE_DNAm", "s_2cell.WE", "s_2cell.WE_DNAm", "s_4cell.WE", "s_Molura.WE", "s_Blast.WE"))
# setnames(sl,3:8, str_replace(colnames(sl)[3:8],'genes.',''))
sl$ex[sl$int] = FALSE
tot.norm = round(1e7/tot.counts,2)
names(tot.norm)[names(tot.norm) == "s_Morula.WE"] <- "s_Molura.WE"


# pdf(file.path(outpath, 'Splicing.dist.pdf'), width=8, height=16)
# 	ggplot(sl, aes(dist)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000) + ggtitle('all')
# 	ggplot(subset(sl, sel.int==TRUE), aes(dist)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('1cell int sort')
# 	ggplot(subset(sl, sel.rat==TRUE), aes(dist)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('1cell int/ex rat sort')
# 	ggplot(subset(sl, sel.samp.rat==TRUE), aes(dist)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('1cell int/ex ~ MII int/ex rat sort')
# 	ggplot(subset(sl, sel.mii==TRUE), aes(dist)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('MII removed')
# 	ggplot(subset(sl, sel.mii.ex==TRUE), aes(dist)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('MII removed ex')
# 	ggplot(subset(sl, sel.mii.int==TRUE), aes(dist)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('MII removed')
# dev.off()
# 
# 
# pdf(file.path(outpath, 'Splicing.dist2.pdf'), width=8, height=16)
# 	ggplot(sl, aes(dist.2)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000) + ggtitle('all')
# 	ggplot(subset(sl, sel.int==TRUE), aes(dist.2)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('1cell int sort')
# 	ggplot(subset(sl, sel.rat==TRUE), aes(dist.2)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('1cell int/ex rat sort')
# 	ggplot(subset(sl, sel.samp.rat==TRUE), aes(dist.2)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('1cell int/ex ~ MII int/ex rat sort')
# 	ggplot(subset(sl, sel.mii==TRUE), aes(dist.2)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('MII removed')
# 	ggplot(subset(sl, sel.mii.ex==TRUE), aes(dist.2)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('MII removed ex')
# 	ggplot(subset(sl, sel.mii.int==TRUE), aes(dist.2)) + geom_histogram(binwidth=10) + facet_grid(samp ~ .) + xlim(-100, 1000)+ ggtitle('MII removed int')
# dev.off()


# pdf(file.path(outpath, 'Splicing.ReadCount.pdf'))

# d = sl[,list(ex=sum(ex), int=sum(int)),by=samp]
# m = melt(d[!str_detect(d$samp,'PA'),])
# m$norm = m$value*tot.norm[m$samp]
# ggplot(m, aes(x=samp, y=norm, fill=variable)) +
# geom_bar(stat="identity", position=position_dodge()) +
# theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
# ggtitle('all')

# d = sl[sl$sel.int, list(ex=sum(ex), int=sum(int)),by=samp]
# m = melt(d[!str_detect(d$samp,'PA'),])
# m$norm = m$value*tot.norm[m$samp]
# ggplot(m, aes(x=samp, y=norm, fill=variable)) +
# geom_bar(stat="identity", position=position_dodge()) +
# theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
# ggtitle('1cell int sort')	

# d = sl[sl$sel.rat, list(ex=sum(ex), int=sum(int)),by=samp]
# m = melt(d[!str_detect(d$samp,'PA'),])
# m$norm = m$value*tot.norm[m$samp]
# ggplot(m, aes(x=samp, y=norm, fill=variable)) +
# geom_bar(stat="identity", position=position_dodge()) +
# theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
# ggtitle('1cell int/ex rat sort')

# d = sl[sl$sel.samp.rat, list(ex=sum(ex), int=sum(int)),by=samp]
# m = melt(d[!str_detect(d$samp,'PA'),])
# m$norm = m$value*tot.norm[m$samp]
# ggplot(m, aes(x=samp, y=norm, fill=variable)) +
# geom_bar(stat="identity", position=position_dodge()) +
# theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
# ggtitle('1cell int/ex ~ MII int/ex rat sort')

d = sl[, list(ex=sum(ex), int=sum(int)),by=samp]
m = melt(d[!str_detect(d$samp,'PA'),])
m$norm = m$value*tot.norm[m$samp]
ggplot(m, aes(x=samp, y=norm, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
  ggtitle('MII removed') +
  ggsave(paste0(outpath, "/Fugaku_no_gene_sel_plot1.pdf"), width = 9.93, height = 5.86)

m.rat = m[,norm[variable=='int']/norm[variable=='ex'], by=samp]
setnames(m.rat,2,'ratio')
ggplot(m.rat, aes(x=samp, y=ratio)) +
  geom_bar(stat="identity") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
  ggtitle('MII removed, Intron/Exon ratio') +
  ggsave(paste0(outpath, "/Fugaku_no_gene_sel_plot2.pdf"), width = 9.93, height = 5.86)


# d = sl[sl$sel.mii.ex, list(ex=sum(ex), int=sum(int)),by=samp]
# m = melt(d[!str_detect(d$samp,'PA'),])
# m$norm = m$value*tot.norm[m$samp]
# ggplot(m, aes(x=samp, y=norm, fill=variable)) +
# geom_bar(stat="identity", position=position_dodge()) +
# theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
# ggtitle('MII removed ex')		

# d = sl[sl$sel.mii.int, list(ex=sum(ex), int=sum(int)),by=samp]
# m = melt(d[!str_detect(d$samp,'PA'),])
# m$norm = m$value*tot.norm[m$samp]
# ggplot(m, aes(x=samp, y=norm, fill=variable)) +
# geom_bar(stat="identity", position=position_dodge()) +
# theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
# ggtitle('MII removed int')		
# dev.off()


# tot.norm = round(1e7/tot.counts,2)
# pdf(file.path(outpath, 'Splicing.ReadCount.MII.removed.pdf'), width=9.93, height=5.86)

# d = sl[sl$sel.mii, list(ex=sum(ex), int=sum(int)),by=samp]
# m = melt(d[!str_detect(d$samp,'PA'),])
# m$norm = m$value*tot.norm[m$samp]
# ggplot(m, aes(x=samp, y=norm, fill=variable)) +
# geom_bar(stat="identity", position=position_dodge()) +
# theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
# ggtitle('MII removed')

m.rat = m[,norm[variable=='int']/norm[variable=='ex'], by=samp]
# setnames(m.rat,2,'int.ex.ratio')
ggplot(m.rat, aes(x=samp, y=V1,fill='int.ex.ratio')) +
  geom_bar(stat="identity", color='black', fill='black') +
  theme_bw() + 
  theme(panel.border = element_rect(size = .1, colour = "black"),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle=90),
        axis.title = element_text(size = 20),
        panel.grid = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position="none") +
  ylab('Unspliced / Spliced') +
  xlab('') +
  ggsave(paste0(outpath, "/Fugaku_no_gene_sel_plot3.pdf"), width = 9.93, height = 5.86)


# dev.off()


sl = data.frame(t(sapply(l, function(x)table(x$don, x$acc))))
colnames(sl) = c('None','don','acc','don-acc')
write.table(sl[,-1], file.path(outpath, 'Splicing_DonorAcceptorAnal.txt'), row.names=T, col.names=T, quote=F, sep='\t')
cat('Thank you, and goodbye!...\n')
