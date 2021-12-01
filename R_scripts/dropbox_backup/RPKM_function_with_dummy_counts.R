counts_RPKM <- function(se, TxDB){
  ebg <- exonsBy(TxDB, by = "gene")
  counts_se <- assay(se)
  counts_se <- counts_se + 1
  lib_size <- colSums(counts_se) 
  ncounts <- t(t(counts_se) / lib_size)
  
  exon_lengths <- width(ebg) 
  exon_lengths_by_gene <- split(exon_lengths, names(exon_lengths))
  gene_lengths <- sapply(exon_lengths_by_gene, sum) 
  names(gene_lengths) <- sub("\\.\\d+", "", names(gene_lengths))
  
  common_genes <- intersect(row.names(ncounts), names(gene_lengths)) 
  subset_ncounts <- ncounts[row.names(ncounts) %in% common_genes, ] 
  gene_lengths <- gene_lengths[names(gene_lengths) %in% common_genes] 
  ncounts <- subset_ncounts/gene_lengths
  rpkm <- ncounts * 1e9 
  rpkm <- as.data.frame(rpkm)
  return(rpkm)  
}
rpkm_new <- counts_RPKM(se, TxDB = TxDb.Mmusculus.UCSC.mm10.knownGene)
rpkm_old <- counts_RPKM(se_PA, TxDB = TxDb.Mmusculus.UCSC.mm9.knownGene)
