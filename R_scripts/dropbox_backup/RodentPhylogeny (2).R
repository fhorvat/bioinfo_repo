### INFO: Given a set of animals and genes it collects the sequences into a fasta format
### DATE: 27.03.2013
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
library(data.table)
library(biomaRt)
library(Biostrings)

#/{1} LIBRARIES


# {2} CODE
	# {{1}} FUNCTIONS
	
	GetIDs = function(id, orgs){
	
		orgs = orgs[names(orgs) != 'human']
		attrs = paste(orgs, '_homolog_ensembl_gene', sep='')
		l = length(attrs)
		set = rep(1:((l %/% 3)+1), each=3)[1:l]

		mart=useMart("ensembl",dataset='hsapiens_gene_ensembl')
		hbm = getBM(attributes=c('ensembl_gene_id'), filters=c('hgnc_symbol'), values=c(id), mart=mart)
		names(hbm) = 'hsapiens'
		l.bm = list()
		for(i in unique(set)){
			print(i)
			l.bm[[i]] = unlist(getBM(attributes=c(attrs[i==set]), filters=c('hgnc_symbol'), values=c(id), mart=mart))
		}
		bm = unlist(l.bm)
		names(bm) = str_replace(names(bm),'_.+','')
		bm = c(hbm, bm)
		return(bm)
	}
	#/{{1}} FUNCTIONS
	
	
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		outpath='/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/MTEvolution/RodentPhylogeny/Sequences/Rodents'
		
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS
		animals=c('guinea pig' = 'cporcellus',
				  'rat' = 'rnorvegicus',
				  'mouse' = 'mmusculus',
				  'rabbit' = 'ocuniculus',
				  # 'squirrle' = 'itridecemlineatus',
				  'kangaroo rat' = 'dordii',
				  'human' = 'hsapiens',
				  'pig' = 'sscrofa'
				  # 'elephant' = 'lafricana',
				  # 'bat' = 'mlucifugus',
				  # 'platypus' = 'oanatinus',
				  # 'dog' = 'cfamiliaris',
				  # 'cow' = 'btaurus'
				  )
		
		hgnc = c('Dmp1') 
				
		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	
	
	l.seq = list()
	### loops over the genes
	for(i in hgnc){
		print(i)
		outdir = file.path(outpath, i)
			dir.create(outdir, showWarnings=F)
		ids = GetIDs(i,animals)
		ids = ids[!duplicated(ids)]
		
		### loops over the organismas and fetches the sequences
		l.data=list()
		for(j in names(ids)){
			
			print(j)
			mart=useMart("ensembl",dataset=paste(j, '_gene_ensembl', sep=''))
			bm = getBM(attributes=c('chromosome_name','start_position','end_position','strand','ensembl_gene_id','coding'), filters=c('ensembl_gene_id'), values=ids[j], mart=mart)
			l.data[[j]] = bm
		}
		d.seq = do.call(rbind, l.data)
		
		### prints the sequences
		if(!nrow(d.seq) == 0){
			
			s = DNAStringSet(d.seq[['coding']])
			names(s) = str_replace(rownames(d.seq),'(_|\\.).+','')
			s = s[order(-width(s))]
			s = s[!duplicated(names(s))]
			names(s) = paste(names(s), i, sep='_')
			
			writeXStringSet(s, file.path(outdir, paste(i,length(s),'fa', sep='.')), format='fasta')
		}
	}

	#/{{3}} MAIN
#/{2} CODE


### takes only 3'rd base

