### INFO: %INFO
### DATE: %DATE
### AUTHOR: %AUTH
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path=file.path(Sys.getenv('HOME'),'bin/MyLib/RFun')
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'ScanLib.R'))
library(data.table)
library(stringr)
library(doMC)

#/{1} LIBRARIES


# {2} CODE
	# {{1}} FUNCTIONS
	#/{{1}} FUNCTIONS
	
	
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS
		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	#/{{3}} MAIN
#/{2} CODE



