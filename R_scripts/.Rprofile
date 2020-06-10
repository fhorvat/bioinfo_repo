### INFO: .Rprofile defining few usefull functions and options
### DATE: Tue Sep 18 14:03:03 2018
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES

######################################################## READ DATA

######################################################## MAIN CODE
### set new hidden enviorment
.env <- new.env()

### adds variable called "pwd" to new enviorment
.env$pwd <- getwd()

### add functions to new enviorment
## sets screen to max. width
.env$wideScreen <- function(howWide = Sys.getenv("COLUMNS")) {
  options(width = as.integer(howWide))
}

## full head of tibble or data.table
.env$headt <- function(df, n = 5) {
  head(as.data.frame(df), n = n)
}

## head and tails
.env$ht <- function(df, n = 5) {
  as.data.frame(rbind(head(df, n = n), tail(df, n = n)))
}

## converts to data.frame
.env$asdf <- function(df) {
  as.data.frame(df)
}

## prints vector in one column
.env$prnt <- function(x) {
  cat(x, sep = "\n")
}


### INFO: parses arguments from command line into named vector
### DATE: 17. 06. 2019.
### AUTHOR: Filip Horvat
.env$parseCommandLineArguments <- function(args) {
  
  # libraries
  require(stringr)
  require(magrittr)
  require(purrr)

  # get named list of arguments
  args_list <-
    args %>%
    unlist(.) %>%
    stringr::str_c(., collapse = " ") %>%
    stringr::str_split(., "--") %>%
    unlist(.) %>%
    .[-1] %>% 
    purrr::map(., function(x){
      
      # split input
      input_split <- 
        x %>% 
        stringr::str_squish(.) %>% 
        stringr::str_split(., " ") %>% 
        unlist(.)
      
      # get arguments
      input_args <- 
        input_split[-1] %>% 
        list(.) %>% 
        set_names(., input_split[1])
      
      # return
      return(input_args)
      
    }) %>% 
    unlist(., recursive = F)
}


### INFO: reads and cleans STAR Log.final.out mapping statistics table
### DATE: 05. 10. 2017.
### AUTHOR: Filip Horvat
.env$read_STAR.Log.final.out <- function(path, reshape = FALSE) {
  
  # libraries
  require(dplyr)
  require(stringr)
  require(magrittr)
  require(readr)
  require(tidyr)

  # read and clean table
  STAR_df <-
    suppressWarnings(readr::read_delim(file = path, delim = "\t", col_names = c("stat", "value"), skip = 5, trim_ws = TRUE)) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(stat = stringr::str_replace(stat, " \\|$", ""),
                  stat = stringr::str_replace(stat, "Number of reads m", "M"),
                  stat = stringr::str_replace(stat, " reads number", ""),
                  value = stringr::str_replace(value, "%", ""),
                  value = as.numeric(value))

  # calculate number of unmapped reads
  columns_mapped <- c("Uniquely mapped",  "Mapped to multiple loci")
  STAR_df %<>%
    dplyr::summarize(value = value[stat == "Number of input reads"] - sum(value[stat %in% columns_mapped]),
                     stat = "Unmapped reads") %>%
    dplyr::bind_rows(STAR_df) %>%
    dplyr::select(stat, value)

  # reshape
  if(reshape){

    columns_out <- c(columns_mapped,  "Unmapped reads", "Number of input reads")
    STAR_df %<>%
      dplyr::filter(stat %in% columns_out) %>%
      dplyr::mutate(stat = factor(stat, levels = columns_out)) %>%
      tidyr::spread(stat, value)  %>%
      dplyr::mutate_all(list(as.integer)) %>%
      dplyr::mutate(log_id = path %>% basename(.) %>% stringr::str_replace(., ".Log.final.out", "")) %>%
      magrittr::set_colnames(colnames(.) %>%  tolower(.) %>% str_replace_all(., " ", "_")) %>%
      dplyr::select(ncol(.), dplyr::everything())

  }

  return(STAR_df)

}


### INFO: transforms gtf to GRanges
### DATE: 05. 04. 2017.
### AUTHOR: vfranke, modified by Filip Horvat
.env$gtfToGRanges <- function(gtf, filter = NULL){

  # libraries
  require(dplyr)
  require(stringr)
  require(magrittr)
  require(tibble)
  require(GenomicRanges)

  gtf %<>%
    magrittr::set_colnames(., c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) %>%
    dplyr::mutate_at(.vars = vars("start", "end"), .funs = list(as.numeric))

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

### INFO: expands GRanges up- and down-stream
### DATE: 2018-05-30
### AUTHOR: Devon Ryan, https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
.env$expandRange <- function(x, upstream=2000, downstream=1000) {
  
  strand_is_minus <- strand(x) == "-"
  on_plus <-- which(!strand_is_minus)
  on_minus <- which(strand_is_minus)
  
  start(x)[on_plus] <- start(x)[on_plus] - upstream
  start(x)[on_minus] <- start(x)[on_minus] - downstream
  
  end(x)[on_plus] <- end(x)[on_plus] + downstream
  end(x)[on_minus] <- end(x)[on_minus] + upstream
  
 return(x)

}

### attaches hidden enviorment
attach(.env)

### fortunes on start
if(interactive())
  try(fortunes::fortune(), silent = TRUE)



