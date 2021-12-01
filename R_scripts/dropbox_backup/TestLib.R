### INFO: Library for testing functions - mostly beta version, or have/can be upgreaded
### DATE: 14.02.2010.
### AUTHOR: v+


# {1} CODE
		#{{ 1 }}
		# maxgap = 0L
		# minoverlap = 1L
		# decreasing=T
		# value='p.value'

		Overlapper = function(tab.tmp){
					
			if(nrow(tab.tmp) != 1){
				tab.tmp = tab.tmp[order(tab.tmp[[value]], decreasing=decreasing),]
				tab.tmp = tab.tmp[!duplicated(tab.tmp$set),]
				pos = c(unique(tab.tmp[,1]), as.character(range(tab.tmp[,2:3])))
				names(pos) = c('chr','start','end')
				val = unlist(tab.tmp[,-c(1:3)])
				str_sub(names(val), -1, -1) = str_c('.',val[str_detect(names(val), 'set')])
				val = val[!str_detect(names(val), 'set')]
				return(quickdf(c(pos,val)))
					
			}else{
				names(tab.tmp)[-c(1:3)] = str_c(names(tab.tmp)[-c(1:3)], tab.tmp$set, sep='.')
				tab.tmp = tab.tmp[,!str_detect(names(tab.tmp), 'set')]
				return(tab.tmp)
			}
		}
		
		MultiSetOverlapper = function(tab, maxgap = 0L, minoverlap = 1L, decreasing=T, value='p.value'){
		
			library(GenomicRanges)
			library(plyr)
			library(stringr)
			tab$chr = as.character(tab$chr)
			
						
			cat('Constructing the ranges...\n')
			ranges = GRanges(
							  seqnames = tab$chr,
							  ranges   = IRanges(start=as.integer(tab$start), end=as.integer(tab$end)),
							)
				
			r = IRanges::reduce(ranges)				 
			cat('Ranges constructed\n\n')
				
				
			cat('Finding overlaps...\n')
			o = as.matrix(findOverlaps(r, ranges, minoverlap=minoverlap, maxgap=maxgap))
			o = o[order(o[,1]),]
			s = split(o, o[,1])
			s = lapply(s, function(x)x[((length(x)/2)+1):length(x)])
			cat('Overlaps found\n\n')
					
			l.data.ovlap = lapply(s, function(i)Overlapper(tab[i,]))
			# l.data.ovlap = foreach(i = 1:length(s)) %dopar% Overlapper(tab[s[[i]],])
			cat('Overlapped sets constructed\n\n')

			# l = lapply(l.data.ovlap, quickdf)
			cat('RbindFill...\n')
			data.ovlap = do.call(plyr::rbind.fill, l.data.ovlap)
			cat('Data frame constructed\n\n')
			
			
			cat('Returning the data...\n')
			return(data.ovlap)
		}
		#{{/1 }}
		
		
#/{1} CODE



