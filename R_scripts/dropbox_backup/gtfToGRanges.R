### INFO: transforms gtf to GRanges
### DATE: 05. 04. 2017.
### AUTHOR: vfranke, modified by Filip Horvat

######################################################## LIBRARIES
library(dplyr)
library(stringr)
library(magrittr)
library(tibble)
library(GenomicRanges)

######################################################## FUNCTIONS
gtfToGRanges <- function(gtf, filter = NULL){
  
  gtf %<>%
    magrittr::set_colnames(., c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) %>%
    dplyr::mutate_at(.vars = vars("start", "end"), .funs = funs(as.numeric(.)))
  
  # gtf has to have 9 columns
  if(ncol(gtf) != 9)
    stop("Number of columns does not match gtf format")
  
  # width of elements have to be positive
  if(any(as.integer(gtf$end) < as.integer(gtf$start))){
    warning("gtf file contains ranges with negative widths...")
    gtf %<>%
      dplyr::filter(end > start)
  }
  
  # filter - for example: exon/transcript/gene/CDS...
  if(!is.null(filter)){
    if(filter %in% gtf$feature){
      gtf %<>%
        dplyr::filter(feature == filter)
    }else{
      stop("The given feature is not present in the gtf file")
    }
  }
  
  ### get exon/transcript/gene ID for the 9th column of gtf table
  # split stings in 9th column, split list again by length (not all genes have all info in that column)
  s <- vector("list", length = length(gtf$attribute))
  s <- stringr::str_extract_all(gtf$attribute, pattern = "gene_id .*?;|transcript_id .*?;|exon_id .*?;|gene_biotype .*?;")
  z <- sapply(s, length)
  a <- split(s, z)
  gtf <- gtf[order(z), ]
  
  # create tables from 9th column with different number of columns in a list
  l <- vector("list", length = length(a))
  l <- lapply(a, function(x){
    
    m <-
      unlist(x, use.names = F) %>% 
      stringr::str_replace_all(., "^.+? |\"|;", "") %>%
      matrix(., ncol = length(x[[1]]), byrow = T) %>%
      magrittr::set_colnames(str_replace(x[[1]], " .+$", "")) %>%
      tibble::as_tibble(.)
    
    return(m)
    
  })
  
  # bind all tables from list in one with filling empty data with NA
  ids <- dplyr::bind_rows(l)
    
  # set all strands which are not +/- to *
  gtf %<>%
    mutate(strand = replace(strand, !(strand %in% c("+", "-")), "*"))
  
  # create and output GRanges
  granges <- GRanges(seqnames = gtf$seqnames,
                     ranges = IRanges(gtf$start, gtf$end),
                     strand = gtf$strand,
                     frame = gtf$frame,
                     feature.type = gtf$feature,
                     .id = 1:nrow(gtf))
  
  values(granges) <- cbind(values(granges), DataFrame(ids[granges$.id, ]))
  values(granges)$.id <- NULL
  
  # return GenomicRanges object
  return(granges)
  
}


