GffToGRanges <- function(gff, filter=NULL){

                library(plyr)
                if(ncol(gff) != 9)
                        stop("Number of columns does not match gff format")

                if(any(gff[,5] < gff[,4])){
                        warning("gff file contains ranges with negative widths...")
                        gff = gff[gff[,5] > gff[,4],]
                }

                if(!is.null(filter)){
                        if(filter %in% gff[,3]){
                                cat("Filtering", filter, "features...\n")
                                gff = gff[gff[,3] == filter,]
                        }else{
                                stop("The given feature is not present in the gff file")
                        }
                }


                cat('Getting the feature ids...\n')
                s = strsplit(gff$V9, split=';')
                z = sapply(s, length)
                a = split(s, z)
                gff = gff[order(z),]
                l = lapply(a, function(x){
                                                d = sub('^ ','', unlist(x, use.names=F))
                                                d = sub('^.+? ','',d)
                                                m = matrix(d, ncol = length(x[[1]]), byrow=T)
                                                colnames(m) = sub(' .+$','',sub('^ ','', x[[1]]))
                                                m})
                ids = rbind.fill(lapply(l, data.frame))

                cat('Constructing the granges...\n')
                gff$V7[!gff$V7 %in% c('+','-')] = '*'
                granges = GRanges(seqnames = gff[,1],
                                                  IRanges(gff[,4],gff[,5]),
                                                  strand = gff[,7],
                                                  frame = gff[,8],
                                                  feature.type = gff[,3],
                                                  .id = 1:nrow(gff))

                values(granges) = cbind(values(granges), DataFrame(ids)[granges$.id,])
                values(granges)$.id = NULL
                return(granges)
        }

