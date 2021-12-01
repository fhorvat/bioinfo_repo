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
library(Cairo)
library(rphast)

#/{1} LIBRARIES


# {2} CODE
	# {{1}} FUNCTIONS
	#/{{1}} FUNCTIONS
	
	
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		
		mammal.tree = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/AccessoryData/UCSCtree/46.mammal.nw'
		
		outpath = '/common/USERS/vfranke/Work/Dicer_RNASeq_29062012/Results/MTEvolution/RodentPhylogeny/NJTree/UCSCTree'
		
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS

		species = c(Hsapiens = 'hg19',
					Mmusculus= 'mm9',
					Rnorvegicus='rn4',
					Dordii='dipOrd1',
					Cporcellus='cavPor3',
					Ocuniculus='oryCun2',
					Btaurus='bosTau4')
		
		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	tree = read.newick.tree(mammal.tree)
	msa.tree = prune.tree(tree, species, all.but=TRUE)
	for(i in names(species))
		msa.tree = str_replace(msa.tree, species[i], i)
	
	treefile=file.path(outpath, 'UCSCtree.nw')
	cat(msa.tree, file=treefile)
	
	rootfile = str_replace(treefile, 'nw', 'Hsap.nw')
	nw_reroot UCSCtree.nw "Hsapiens" > UCSCtree.Hsap.nw
	nw_display -s -w 600 -c css.map -b 'opacity:0' -l 'font-size:large;font-family:bold' -d 'stroke-width:4' 'UCSCtree.Hsap.nw' > 'UCSCtree.Hsap.svg'
	
	#/{{3}} MAIN
#/{2} CODE


