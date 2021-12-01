### INFO: 
### DATE: Mon May 14 01:18:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/test")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# rmsk.txt.gz path
rmsk.txt.gz_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/rmskVIZ.mm10.20180512.raw.txt.gz"

# rmsk.fa.out.gz path
rmsk.fa.out.gz_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/rmsk.mm10.20181405.raw.fa.out.gz"

######################################################## READ DATA
# bin	      607	      smallint(5) unsigned	Indexing field to speed chromosome range queries.
# swScore	  12937	    int(10) unsigned	Smith Waterman alignment score
# milliDiv	106	      int(10) unsigned	Base mismatches in parts per thousand
# milliDel	13	      int(10) unsigned	Bases deleted in parts per thousand
# milliIns	2	        int(10) unsigned	Bases inserted in parts per thousand
# genoName	chr1	    varchar(255)	Genomic sequence name
# genoStart	3000000	  int(10) unsigned	Start in genomic sequence
# genoEnd	  3000097	  int(10) unsigned	End in genomic sequence
# genoLeft	192471874	int(11)	-#bases after match in genomic sequence
# strand	  -	        char(1)	Relative orientation + or -
# repName	  L1_Mus3	  varchar(255)	Name of repeat
# repClass	LINE	    varchar(255)	Class of repeat
# repFamily	L1	      varchar(255)	Family of repeat
# repStart	3487	    int(11)	Start (if strand is +) or -#bases after match (if strand is -) in repeat sequence
# repEnd	  3592	    int(11)	End in repeat sequence
# repLeft	  3055	    int(11)	-#bases after match (if strand is +) or start (if strand is -) in repeat sequence
# id	      1	        int(11)	First digit of id field in RepeatMasker .out file. Best ignored.

# read rmsk.txt.gz
rmsk_txt <- readr::read_table2(file = rmsk.txt.gz_path, col_names = c("bin", "swScore", "milliDiv", "milliDel", "milliIns", 
                                                                      "genoName", "genoStart", "genoEnd", "genoLeft", "strand", 
                                                                      "repName", "repClass", "repFamily", "repStart", "repEnd", 
                                                                      "repLeft", "id")) %>% 
  dplyr::mutate(strand = replace(strand, strand == "C", "-"), 
                genoStart = genoStart + 1)

# read rmsk.out
rmsk_out <- 
  readr::read_table2(file = rmsk.fa.out.gz_path, skip = 3, col_names = c("swScore", "percDiv", "percDel", "percIns", 
                                                                         "genoName", "genoStart", "genoEnd", "genoLeft", "strand", 
                                                                         "repName", "repClass_repFamily", "repStart", "repEnd", 
                                                                         "repLeft", "id")) %>% 
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-"))


######################################################## MAIN CODE
# filter and arrange
rmsk_out_filt <- 
  rmsk_out %>% 
  dplyr::select(genoName, genoStart, genoEnd, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id) %>% 
  dplyr::mutate(repStart = str_remove_all(repStart, "\\(|\\)") %>% as.integer(.), 
                repLeft = str_remove_all(repLeft, "\\(|\\)") %>% as.integer(.)) %>% 
  dplyr::arrange(genoName, genoStart) 

# repEnd can be thought of as the "middle" coordinate of the three.  It is always 
# a positive integer regardless of the strand the repeat element hits.  It 
# represents the end coordinate of the matching part of the repeat element 
# in "repeat element coordinates".   The "beginning" coordinate of the 
# matching part of a repeat element is repStart (for + oriented hits) or 
# repLeft (for - oriented hits). 

# For + oriented hits, repStart & repEnd are the proper coordinates of the 
# matching part of the element.  RepLeft in this case is a numerical value 
# which may be used (via the equation repEnd-repLeft) to obtain the size 
# of the repeat element ("Left" in the sense of "repeat remaining" unaligned). 
# 
# For - oriented hits, repLeft & repEnd are the proper coordinates of the 
# matching part of the element.  RepLeft will be upstream from (smaller than) 
# repEnd for these neg-strand alignments.  RepStart in this case is 
# a numerical value which may be used (via the equation repEnd-repStart) 
# to obtain the size of the repeat element. 

rmsk_line <-
  rmsk_out_filt %>%
  dplyr::filter(repClass == "LINE", repFamily == "L1")

line_plus <-
  rmsk_line %>%
  dplyr::filter(strand == "+", repStart < 50, repLeft < 50) %>%
  dplyr::filter(repName == "L1Mb_T")

line_minus <-
  rmsk_line %>%
  dplyr::filter(strand == "-", repLeft < 50, repStart < 50)
