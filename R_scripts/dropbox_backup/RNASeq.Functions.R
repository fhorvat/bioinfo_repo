# INFO: function for analysis of rnaseq data
# DATE: 20.05.2011.
# AUTH: v+

# {1} FUNCTIONS
	
	# {{1}}
	# calculates rpkm for a given vector of counts and transcripts
	RPKM = function(expr, width, readnum){
	
		rpkm = 10^9 * expr /(width * readnum)
		return(rpkm)
	}
	#/{{1}}
	
	
	# {{2}}
	# does DESeq differential expression on a table of counts
	# for no replicates
	DiffExpNoRep = function(expr){
	
		library(DESeq)
		cds = newCountDataSet(expr, c('T','N'))
		cds = estimateSizeFactors(cds)
		cds = estimateVarianceFunctions(cds, method="blind" )
		d = nbinomTest(cds, 'T', 'N')
		d$padj[is.na(d$padj)] = 1
		return(d)
	}
	
	# one sample has replicates
	DiffExpOneRep = function(expr){
	
		library(DESeq)
		cds = newCountDataSet(expr, c('T','N'))
		cds = estimateSizeFactors(cds)
		cds = estimateVarianceFunctions(cds, method="blind" )
		d = nbinomTest(cds, 'T', 'N')
		d$padj[is.na(d$padj)] = 1
		return(d)
	}
	#/{{2}}
	
#/{1} FUNCTIONS
