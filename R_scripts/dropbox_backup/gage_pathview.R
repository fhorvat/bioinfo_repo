gagePathview <- function(x, y){
  # gage
  deseq2.fc <- res$log2FoldChange
  names(deseq2.fc) <- rownames(res)
  fc.kegg.p <- gage(deseq2.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
  
  # upregulated pathways
  sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
  path.ids <- rownames(fc.kegg.p$greater)[sel]
  path.ids <- substr(path.ids, 1, 8)
  df_upregulated <- data.frame(q.value = fc.kegg.p$greater[sel, "q.val"])
  if (nrow(df_upregulated) > 1){
    df_upregulated$number <- 1 : nrow(df_upregulated)  
  }
  write.csv(df_upregulated, paste0(x, "_vs_", y, "_", sampleTable_exp$Time.Course[1], "_upregulated_all.csv"))
  pv_upregulated <- sapply(path.ids, function(pid) pathview(gene.data = deseq2.fc, 
                                                            pathway.id = pid, 
                                                            species = "mmu",
                                                            limit = list(gene = 3, cpd = 2),
                                                            bins = list(gene = 20, cpd = 20),
                                                            out.suffix = paste0(x, "_vs_", y,  "_", sampleTable_exp$Time.Course[1], "_upregulated_all")))
  
  # downregulated pathways
  sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
  path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
  path.ids.l <- substr(path.ids.l, 1, 8)
  df_downregulated <- data.frame(q.value = fc.kegg.p$less[sel.l, "q.val"])
  if (nrow(df_downregulated) > 1){
    df_downregulated$number <- 1 : nrow(df_downregulated)  
  }
  write.csv(df_downregulated, paste0(x, "_vs_", y,  "_", sampleTable_exp$Time.Course[1], "_downregulated_all.csv"))
  pv_downregulated <- sapply(path.ids.l, function(pid) pathview(gene.data = deseq2.fc, 
                                                                pathway.id = pid, 
                                                                species = "mmu",
                                                                limit = list(gene = 3, cpd = 2),
                                                                bins = list(gene = 20, cpd = 20),
                                                                out.suffix = paste0(x, "_vs_", y, "_", sampleTable_exp$Time.Course[1], "_downregulated_all")))
  print("end pathview output")
}