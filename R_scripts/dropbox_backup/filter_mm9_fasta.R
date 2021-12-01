library(BSgenome.Mmusculus.UCSC.mm9)
library(GenomicRanges)
library(rtracklayer)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/references/filtered_fasta/1000bp_overhang")
reads <- import.bed(con = "../../original_genome/L1_all_20160418.bed.gz")
reads_filtered <- reads[width(reads) > 6000]

# fasta
L1_fasta <- getSeq(Mmusculus,
                   names = as.character(seqnames(reads_filtered)), 
                   start = start(reads_filtered) - 1000, 
                   end = end(reads_filtered) + 1000, 
                   strand = as.character(strand(reads_filtered)))
names(L1_fasta) <- paste0(as.character(seqnames(reads_filtered)), 
                          ":", 
                          start(reads_filtered) - 1000,
                          "-", 
                          end(reads_filtered) + 1000,
                          ":",
                          as.character(strand(reads_filtered)))
writeXStringSet(L1_fasta, "L1_6000bp_1000bp_overhang.fasta")

# 6-column bed from fasta
L1_bed <- data.frame(chrom = names(L1_fasta),
                     chromStart = 0, 
                     chromEnd = width(L1_fasta),
                     name = mcols(reads_filtered)$name,
                     score = 500,
                     strand = as.character(strand(reads_filtered)))

write.table(L1_bed, "L1_6000bp_1000bp_overhang_from_fasta.bed", quote = F, sep = "\t", row.names = F, col.names = F)

