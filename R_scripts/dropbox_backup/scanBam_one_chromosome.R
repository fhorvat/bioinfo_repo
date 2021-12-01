lib_path <- "/common/WORK/fhorvat/code_library/R_scripts/vfranke"
source(file.path(lib_path, "BamWorkers.R"))

bam_file <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/common/GarciaLopez_2014_development_smallRNA_GSE45983/Data/Mapped/STAR_mm10/s_ooctye_Aligned.sortedByCoord.out.bam"

# find chromosomes 
chrs <- chrFinder(bam_file)
chrs <- chrs[chrs$chr != "chrM",]
chr <- "chr1"
which_ranges <- GRanges(chr, IRanges(1, chrs$chrlen[chrs$chr == chr]))

param <- Rsamtools::ScanBamParam(which = which_ranges, what = c("rname", "strand", "pos", "seq", "cigar"))
bam_in <- Rsamtools::scanBam(bam_file, param = param)
