library("DECIPHER")
library("Biostrings")

retrotransposon <- readDNAStringSet("Sampled_MTs_all_classes_20160801-223646.fasta", "fasta")
retrotransposon_names <- unlist(lapply(X = strsplit(names(retrotransposon), split = "\\|"), FUN = function(X) X[2]))
retrotransposon_list <- split(retrotransposon, rep(1:20, each = 200))
retrotransposon_list_oriented <- lapply(X = retrotransposon_list, 
                                        FUN = function(X){
                                          OrientNucleotides(myXStringSet = X, 
                                                            reference = which.max(width(X)), 
                                                            type = "sequences", 
                                                            orientation = "all", 
                                                            threshold = 0.05, 
                                                            verbose = TRUE, 
                                                            processors = 4)
                                        })

names(retrotransposon_list_oriented) <- NULL 
retrotransposon_oriented <- do.call(c, retrotransposon_list_oriented)
writeXStringSet(retrotransposon_oriented, "Sampled_MTs_all_classes_20160801-223646_oriented.fasta")