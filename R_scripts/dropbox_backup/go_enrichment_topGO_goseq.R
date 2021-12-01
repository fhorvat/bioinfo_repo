library("dynamicTreeCut")
library("WGCNA")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("topGO")
library("plyr")
library("goseq")


results_all_no_NA <- results_all[rowSums(is.na(results_all)) < 1, ]
results_all_output_distance <- as.matrix(results_all_no_NA)
results_all_output_distance <- dist(results_all_output_distance)
hc.rows <- hclust(results_all_output_distance, method = "ward.D")

# Dynamic Tree Cut
maxCoreScatter <- 0.93
minGap <- (1 - maxCoreScatter) * 3/4
dynamicCut <- cutreeDynamic(hc.rows,
                            method = "hybrid", 
                            distM = as.matrix(results_all_output_distance), 
                            deepSplit = 2, 
                            maxCoreScatter = maxCoreScatter, 
                            minGap = minGap, 
                            maxAbsCoreScatter = NULL, 
                            minAbsGap = NULL)

results_clusters <- data.frame(results_all_no_NA, cluster = dynamicCut)

# GO enrichment per clusters (cluster vs. background)
results_clusters$ensembl_gene_id <- mapIds(org.Mm.eg.db,
                                           keys = row.names(results_clusters),
                                           column ="ENSEMBL",
                                           keytype = "ENTREZID",
                                           multiVals = "first")
results_clusters <- results_clusters[-which(is.na(results_clusters$ensembl_gene_id)), ]
results_clusters <- results_clusters[-which(duplicated(results_clusters$ensembl_gene_id)), ]
rownames(results_clusters) <- results_clusters$ensembl_gene_id
results_clusters <- data.frame(results_clusters$cluster, row.names = rownames(results_clusters))
colnames(results_clusters) <- "cluster"

# list of data.frames (each cluster has its own data.frame)
results_clusters_list <- lapply(seq_len(max(results_clusters$cluster)), function(X) results_clusters)

# adding new ID column to each data.frame with 0 or 1 (one cluster per data.frame get's ID 1)
for (i in 1 : max(results_clusters$cluster)){
  results_clusters_list[[i]]$ID <-  as.numeric(results_clusters_list[[i]]$cluster == i)  
}


# topGO 
topGO_write <- function(results_clusters_list){
  results_vector <- results_clusters_list$ID
  names(results_vector) <- rownames(results_clusters_list)
  results_vector <- factor(results_vector)
  onts <- c("MF", "BP", "CC")
  tab <- as.list(onts)
  names(tab) <- onts
  
  for(j in 1:3){
    
    ## prepare data
    tgd <- new("topGOdata", ontology = onts[j], allGenes = results_vector, nodeSize = 5,
               annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")
    
    ## run tests
    resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")
    
    ## look at results
    tab[[j]] <- GenTable(tgd, Fisher.elim = resultTopGO.elim, 
                         Fisher.classic = resultTopGO.classic,
                         orderBy = "Fisher.classic", topNodes = 20)
    
  }
  topGOResults <- rbind.fill(tab)
  return(topGOResults)
}


topGO_results <- lapply(results_clusters_list, topGO_write)
lapply(1:length(topGO_results), function(i) write.csv(topGO_results[[i]], 
                                                      file = paste0("TopGOResults_cluster", i, "all_genes.csv"),
                                                      row.names = FALSE))

# goseq 
goseq_function <- function(results_clusters_list){
  results_vector <- results_clusters_list$ID
  names(results_vector) <- rownames(results_clusters_list)
  pwf <- nullp(results_vector, "mm10", "ensGene")
  GO.wall <- goseq(pwf, "mm10","ensGene")
  enriched.GO <- GO.wall[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05, ]
  enriched.GO <- enriched.GO[, c(1, 6, 7, 2 : 5)]
  return(enriched.GO)
}

goseq_results <- lapply(results_clusters_list, goseq_function)
lapply(1:length(goseq_results), function(i) write.csv(goseq_results[[i]], 
                                                      file = paste0("goseqResults_cluster", i, "all_genes.csv"),
                                                      row.names = FALSE))

