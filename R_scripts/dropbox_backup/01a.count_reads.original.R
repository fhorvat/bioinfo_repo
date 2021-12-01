setwd("/common/WORK/rosa/MT_elements/Scripts/")

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)

# get command arguments
# read in the name of the data file, name of the outfile, fold of crossvalidaton
args <- commandArgs(TRUE)

workdir=args[length(args)-3]
bamFile=args[length(args)-2]
rangeFile=args[length(args)-1]
outFile=args[length(args)]


print(args)

setwd(workdir)

# read in info for MT elements
mt = read.delim(rangeFile, stringsAsFactors = FALSE)

# make ranges for MT elements
mtRanges = GRanges(IRanges(mt$genoStart, mt$genoEnd), seqnames = Rle(mt$genoName), strand = mt$strand, name = mt$repName)

# check if the experiment is paired end or single end
paired = testPairedEndBam(bamFile)

if (paired){
  reads = readGAlignmentPairs(bamFile, use.names = TRUE)
} else {
  reads = readGAlignments(bamFile, use.names = T)
}

# find all overlaps between reads and ranges
my.readOverlaps = findOverlaps(mtRanges, reads, ignore.strand = T)


my.readsTable = data.frame(ranges = as.character(mtRanges$name[queryHits(my.readOverlaps)]), read_name = names(reads)[subjectHits(my.readOverlaps)])
rm(reads)
gc()

my.uniqueReadsTable =  unique(my.readsTable)

counts =  my.uniqueReadsTable[!(my.uniqueReadsTable$read_name %in% names(which(table(my.uniqueReadsTable$read_name)>1))),]


#save(list = c("reads", "my.readsTable", "my.uniqueReadsTable","counts"), file=outFile)

#outFile2 = paste(outFile, "v2", sep = ".")


final.counts = as.data.frame(table(counts$ranges))

rownames(final.counts) = final.counts[,1]
column = gsub("\\.bam", "", gsub("^.*\\/","",bamFile))
final.counts[,1] <- NULL
colnames(final.counts) = column


# save the results
save(final.counts, file=outFile)
