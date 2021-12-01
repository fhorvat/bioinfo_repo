### INFO: Selects the sequences that overlap each domain
### DATE: 28.03.2013
### AUTHOR: Vedran Franke
rm(list=ls())


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path='/home/members/vfranke/MyLib'
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'ScanLib.R'))
source(file.path(lib.path, 'ExternalApps.R'))
library(data.table)
library(ape)

#/{1} LIBRARIES


# {2} CODE
	# {{1}} FUNCTIONS
	#/{{1}} FUNCTIONS
	
	
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		inpath='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/MTEvolution/RodentPhylogeny/'
		
		seq.path = file.path(inpath, 'Sequences/Rodents')
		mult.path = file.path(inpath, 'MultipleAlignment/ClustalW')
		tree.path = file.path(inpath, 'NJTree')
		
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS
		domains = IRanges(start=c(48, 434, 620, 876, 1283, 1665), end=c(199, 543, 712, 998, 1356, 1830))
		names(domains) = c('DEXDc','HELOCc','dsRNA_bind','PAZ','Rnc','RIBOc')
		
		
		data(BLOSUM62)

		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	### draws the nj tree from the mals
	infiles = list.files(seq.path, full.names=T, pattern='fa$')
	infiles = infiles[Reduce('|', lapply(c('AA','DNA'), function(x)str_detect(infiles,x)))]
	infiles = infiles[-2]
	
	for(i in 1:length(infiles)){
		
		file = infiles[i]
		name = str_replace(basename(file),'.fa','')
		print(name)
		what = ifelse(str_detect(name, 'AA'),'AA','DNA')
		outfile=file.path(mult.path, paste(name, 'ClustalO', 'fa', sep='.'))
		a = ClustalO(file, outfile, what=what, format='fa', force=T)
		names(a) = substr(names(a), 1, 10)
		writeXStringSet(a, outfile, format='fasta')
		
		if(what == 'AA')
			ind = which(apply(as.matrix(a),2, function(x)length(unique(x))==1))

		for(j in 1:(length(domains))){
		
			domain = domains[j]
			dom.name = names(domain)
			print(dom.name)
				
			if(what=='AA'){
				en = endoapply(a, function(x)x[start(domain):end(domain)])
			}else if(str_detect(outfile, 'neut')) {
				ind.t = ind[ind >= start(domain) & ind <= end(domain)]
				en = endoapply(a, function(x)x[ind.t])
			}
			inf= sum(apply(as.matrix(en),2,function(x)length(unique(x)))!=1)
			writeXStringSet(en, file.path(mult.path, paste(dom.name, name,'w',unique(width(en)),'i',inf,'fa', sep='.')), format='fasta')
		}
	}
	#/{{3}} MAIN
#/{2} CODE


