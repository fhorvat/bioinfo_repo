### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression/GO_terms_enrichment")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(biomaRt)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(clusterProfiler)
library(GOplot)
library(org.Mm.eg.db)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# modified GOplot::GOBubble which uses ggrepel::geom_text_repel instead of ggplot2::geom_text
GOBubble2 <- function(data, display, title, colour, labels, ID = T, table.legend = T, table.col = T, bg.col = F, pvalue_threshold = 0.01){
  
  zscore <- adj_pval <- category <- count <- id <- term <- NULL
  
  if(missing(display)){
    display <- "single"
  }
  
  if (missing(title)){
    title <- ""
  }
  
  if(missing(colour)){
    cols <- c("chartreuse4", "brown2", "cornflowerblue")
  }else{
    cols <- colour
  }
  
  if(missing(labels)){
    labels <- 5
  }
  
  if(bg.col == T & display == "single"){
    cat("Parameter bg.col will be ignored. To use the parameter change display to 'multiple'")
  }
  
  colnames(data) <- tolower(colnames(data))
  
  if(!"count" %in% colnames(data)) {
    rang <- c(5, 5)
    data$count <- rep(1, dim(data)[1])
  }else{
    rang <- c(1, 30)
  }
  
  data$adj_pval <- -log(data$adj_pval, 10)
  sub <- data[!duplicated(data$term), ]
  pvalue_threshold_log <- -log(pvalue_threshold, 10)
  
  g <- 
    ggplot(sub, aes(zscore, adj_pval, fill = category, size = count)) +
    labs(title = title, x = "z-score", y = "-log (adj p-value)") +
    geom_point(shape = 21, col = "black", alpha = 1/2) +
    geom_hline(yintercept = pvalue_threshold_log, col = "orange") + 
    scale_x_continuous(limits = c(min(sub$zscore) - 0.2, max(sub$zscore) + 0.2)) +
    scale_size(range = rang, guide = "none")
  
  if(!is.character(labels)){
    sub2 <- subset(sub, subset = sub$adj_pval >= labels)
  }else{
    sub2 <- subset(sub, sub$id %in% labels | sub$term %in% labels)
  }
  
  if (display == "single") {
    
    g <- 
      g + 
      scale_fill_manual("Category", values = c(BP = cols[1], CC = cols[2], MF = cols[3]), labels = c("Biological Process", "Cellular Component", "Molecular Function")) + 
      theme(legend.position = "bottom") +
      annotate("text", x = min(sub$zscore) + 0.2, y = pvalue_threshold_log + 0.1, label = paste0(pvalue_threshold, " p-value threshold"), colour = "orange", size = 4)
    
    if(ID){
      g <- 
        g + 
        ggrepel::geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval, label = id), size = 5)
    } else{
      g <- 
        g + 
        ggrepel::geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 4)
    } 
    
    if((table.legend) & (nrow(sub2) > 0)){
      
      if(table.col){
        table <- GOplot:::draw_table(sub2, col = cols)
      }else{
        table <- GOplot:::draw_table(sub2)
      }
      
      g <- 
        g + theme(axis.text = element_text(size = 14),
                  axis.line = element_line(colour = "grey80"),
                  axis.ticks = element_line(colour = "grey80"),
                  axis.title = element_text(size = 14, face = "bold"),
                  panel.background = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  plot.background = element_blank())
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      grid.arrange(g, table, ncol = 2)
      
    }else{
      
      g <- 
        g + theme(axis.text = element_text(size = 14), 
                  axis.line = element_line(colour = "grey80"),
                  axis.ticks = element_line(colour = "grey80"),
                  axis.title = element_text(size = 14, face = "bold"),
                  panel.background = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  plot.background = element_blank())
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      grid.arrange(g, ncol = 1)
      
    }
    
  }else{
    
    if(bg.col){
      dummy_col <- data.frame(category = c("BP", "CC", "MF"), 
                              adj_pval = sub$adj_pval[1:3], 
                              zscore = sub$zscore[1:3],
                              size = 1:3, 
                              count = 1:3)
      
      g <- 
        g + 
        geom_rect(data = dummy_col, aes(fill = category), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) + 
        facet_grid(. ~ category, space = "free_x", scales = "free_x") + 
        scale_fill_manual(values = c(BP = cols[1], CC = cols[2], MF = cols[3]), guide = "none", drop = F)
      
    }else{
      
      g <- 
        g +
        facet_grid(. ~ category, space = "free_x", scales = "free_x") + 
        scale_fill_manual(values = c(BP = cols[1], CC = cols[2], MF = cols[3]), guide = "none")
    }
    
    if(ID){
      
      g <- 
        g + 
        ggrepel::geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval, label = id), size = 5) + 
        theme(axis.title = element_text(size = 14, face = "bold"), 
              axis.text = element_text(size = 14), 
              axis.line = element_line(colour = "grey80"), 
              axis.ticks = element_line(colour = "grey80"),
              panel.border = element_rect(fill = "transparent", colour = "grey80"), 
              panel.background = element_blank(),
              panel.grid = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_blank())
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      grid.arrange(g, ncol = 1)
      
    }else{
      
      g <- 
        g + 
        ggrepel::geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 5) + 
        theme(axis.title = element_text(size = 14, face = "bold"), 
              axis.text = element_text(size = 14),
              axis.line = element_line(colour = "grey80"),
              axis.ticks = element_line(colour = "grey80"),
              panel.border = element_rect(fill = "transparent", colour = "grey80"),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_blank())
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      grid.arrange(g, ncol = 1)
      
    }
  }
}

######################################################## PATH VARIABLES
### experiment
# set experiment name
experiment <- "Lnc1_KO"
experiment_name <- "Lnc1_KO"

### working dir
# set working directory 
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/diffExp/lnc1_KO"))


### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### other experiment paths
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec")

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- inpath


### documentation
# set ensembl version
ensembl_version <- 93

# sample table path
sample_table_path <- file.path(documentation_path, "lnc1_KO.RNAseq.20181211.sampleTable.clean.csv")

# stats and tracks path
stats_and_tracks_path <- list.files(path = mapped_path, pattern = ".*\\.stats_and_tracks\\.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.se\\.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.FPKM_mean\\.csv$"), full.names = T)

# differential expression results path
diffExp_path <- file.path(analysis_path, str_c("diffExp", "ensembl", ensembl_version,
                                               "lnc1_Null_vs_WT", "DESeq2.protein_coding.significant_results.xlsx", sep = "."))


### genome
# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# GO terms info
go_terms_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.GOterms.csv$"), full.names = T)


######################################################## READ DATA
# read gene info
genes_info <- 
  readr::read_csv(genes_info_path) %>% 
  dplyr::filter(gene_biotype == "protein_coding")

# read GO terms info
go_term_info <- readr::read_csv(go_terms_info_path)

### read list of diff. expressed genes
# set stages
stages <- c("GV", "MII")

# read different stages
diffExp_results <- purrr::map(stages, function(stage){
  
  # read sheet
  read.xlsx(xlsxFile = diffExp_path, sheet = str_c(stage, ".lnc1_Null_vs_WT")) %>% 
    as_tibble(.)
  
}) %>% 
  set_names(stages)

######################################################## MAIN CODE
### do the enrichment for both stages
go_enrich_list <- purrr::map(stages, function(stage){
  
  # get results table for one stage
  diff_df <- diffExp_results[[stage]]
  
  # prepare gene list
  geneList <- diff_df$log2FoldChange
  names(geneList) <- diff_df$gene_id
  geneList <- sort(geneList, decreasing = TRUE)
  
  ### GO terms enrichment
  # do GO enrichment across 3 GO categories
  go_categories <- c("CC", "MF", "BP")
  
  # enrich
  go_enrich <- purrr::map(go_categories, function(go_cat){
    
    # enrich GO terms
    ego <- clusterProfiler::enrichGO(gene = diff_df$gene_id,
                                     universe = genes_info$gene_id,
                                     keyType = "ENSEMBL",
                                     OrgDb = org.Mm.eg.db,
                                     ont = go_cat,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.01,
                                     qvalueCutoff = 0.05,
                                     readable = TRUE)
    
    if(nrow(ego) > 0){
      
      # save to pdf
      pdf(file = file.path(outpath, str_c("GOEnrich.visualization", stage, go_cat, "full.pdf", sep = ".")),
          width = 15,
          height = 10)
      tryCatch({print(barplot(ego, drop = TRUE, showCategory = 15))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      tryCatch({print(dotplot(ego))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      tryCatch({print(emapplot(ego))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      tryCatch({print(goplot(ego, geom = "label"))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      tryCatch({print(cnetplot(ego, categorySize = "pvalue", foldChange = geneList))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      dev.off()
      
      # simplify enriched GO terms by removing redundancy
      ego_simple <- clusterProfiler::simplify(x = ego,
                                              cutoff = 0.7,
                                              by = "p.adjust",
                                              select_fun = min)
      
      # save visualizations
      if(nrow(ego_simple) > 0){
        
        # save to pdf
        pdf(file = file.path(outpath, str_c("GOEnrich.visualization", stage, go_cat, "simple.pdf", sep = ".")),
            width = 15,
            height = 10)
        tryCatch({print(barplot(ego_simple, drop = TRUE, showCategory = 15))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
        tryCatch({print(dotplot(ego_simple))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
        tryCatch({print(emapplot(ego_simple))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
        tryCatch({print(goplot(ego_simple, geom = "label"))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
        tryCatch({print(cnetplot(ego_simple, categorySize = "pvalue", foldChange = geneList))}, error = function(e){cat("ERROR :" , conditionMessage(e), "\n")})
        dev.off()
        
        # convert to tibble
        ego_df_simple <- 
          ego_simple %>% 
          as.data.frame(.) %>% 
          as_tibble(.) %>% 
          dplyr::mutate(GO_category = go_cat) %>% 
          dplyr::select(ID, GO_category, everything())
        
      }
      
      # convert to tibble
      ego_df <- 
        ego %>% 
        as.data.frame(.) %>% 
        as_tibble(.) %>% 
        dplyr::mutate(GO_category = go_cat) %>% 
        dplyr::select(ID, GO_category, everything())
      
    }else{
      
      # set as empty tibble
      ego_df <- tibble()
      ego_df_simple <- tibble()
      
    }
    
    # return
    return(list(full_ego = ego_df, simple_ego = ego_df_simple))
    
  })
  
}) %>% 
  set_names(., stages)

# save results of enrichment
saveRDS(go_enrich_list, file = file.path(outpath, str_c("diffExp", 
                                                        "ensembl", ensembl_version,
                                                        "lnc1_Null_vs_WT",
                                                        "DESeq2.protein_coding.significant_results.GO_terms.RDS", sep = ".")))



### write results, plot bubble plot
purrr::map(stages, function(stage){
  
  # get enrichment table for one stage
  go_enrich <- go_enrich_list[[stage]]
  
  # join results - full GO enrich
  ego_full_df <- 
    purrr::map(go_enrich, function(x) x$full_ego) %>% 
    dplyr::bind_rows(.) %>% 
    dplyr::filter(!duplicated(ID)) %>% 
    dplyr::rename(GO_id = ID)  %T>% 
    readr::write_csv(., path = file.path(outpath, str_c("GOEnrich.results", stage, "full.csv", sep = "."))) %>%
    dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
    dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                  Genes = str_replace_all(Genes, "/", ", "))
  
  # simplified GO enrich
  ego_simple_df <- 
    purrr::map(go_enrich, function(x) x$simple_ego) %>% 
    dplyr::bind_rows(.) %>% 
    dplyr::filter(!duplicated(ID)) %>% 
    dplyr::rename(GO_id = ID) %T>%
    readr::write_csv(., path = file.path(outpath, str_c("GOEnrich.results", stage, "simple.csv", sep = "."))) %>% 
    dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
    dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                  Genes = str_replace_all(Genes, "/", ", "))
  
  # # rename columns in results
  # genes_info_df <-
  #   diffExp_results[[stage]] %>% 
  #   dplyr::left_join(., go_term_info, by = "gene_id") %>% 
  #   dplyr::rename(padj = padj, ID = gene_name, logFC = log2FoldChange)
  # 
  # ### plot results
  # for(result in c("full", "simple")){
  #   
  #   result <- "full"
  #   
  #   # set results table
  #   if(result == "full"){
  #     circ <- GOplot::circle_dat(terms = ego_full_df, genes = genes_info_df)
  #   }else{
  #     circ <- GOplot::circle_dat(terms = ego_simple_df, genes = genes_info_df)
  #   }
  #   
  #   # GO bubble plot with table
  #   png(file = file.path(outpath, str_c("GOEnrich.bubble.single.table", stage, result, "png", sep = ".")),
  #       width = 1800,
  #       height = 1000,
  #       units = "px",
  #       type = "cairo")
  #   GOBubble2(data = circ,
  #             labels = 1,
  #             display = "single",
  #             table.legend = T,
  #             pvalue_threshold = 0.01,
  #             title = str_c("Enriched GO terms ", "(", result, ") - ", stage))
  #   dev.off()
  #   
  #   # GO bubble plot facet
  #   png(file = file.path(outpath, str_c("GOEnrich.bubble.facet", stage, result, "png", sep = ".")),
  #       width = 1000,
  #       height = 1000,
  #       units = "px",
  #       type = "cairo")
  #   GOBubble2(data = circ,
  #             labels = 1,
  #             display = "multiple",
  #             bg.col = T,
  #             pvalue_threshold = 0.01,
  #             title = str_c("Enriched GO terms ", "(", result, ") - ", stage))
  #   dev.off()
  # 
  # }
  
})


### pathview
# KEGG set for mouse
kg.mouse <- kegg.gsets("mouse")
kegg.gs <- kg.mouse$kg.sets[kg.mouse$sigmet.idx]

### write results, plot bubble plot
purrr::map(stages, function(stage){
  
  # get results table for one stage
  diff_df <- diffExp_results[[stage]]
  
  # prepare gene list
  geneList <- diff_df$log2FoldChange
  names(geneList) <- diff_df$gene_id
  geneList <- sort(geneList, decreasing = TRUE)
  
  # gage
  fc.kegg.p <- gage(geneList, gsets = kegg.gs, ref = NULL, samp = NULL)
  
  ### get upregulated pathways
  # upregulated pathways
  sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
  path.ids <- rownames(fc.kegg.p$greater)[sel]
  path.ids <- substr(path.ids, 1, 8)
  
  # create table
  df_upregulated <- 
    tibble(q.value = fc.kegg.p$greater[sel, "q.val"]) %>% 
    dplyr::mutate(number = 1:nrow(.)) %T>% 
    readr::

  # save tab;e
  write.csv(df_upregulated, paste0("./gage_pathview/", "gage_", sample_list, "_upreg.csv"))
  
  
  pv_upregulated <- sapply(path.ids, function(pid) pathview(gene.data = deseq2.fc,
                                                            pathway.id = pid,
                                                            species = "mmu",
                                                            limit = list(gene = 3, cpd = 2),
                                                            bins = list(gene = 20, cpd = 20),
                                                            out.suffix = paste0(sample_list, "_upreg"),
                                                            kegg.dir = "./gage_pathview"))
  
  # downregulated pathways
  sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
  path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
  path.ids.l <- substr(path.ids.l, 1, 8)
  df_downregulated <- data.frame(q.value = fc.kegg.p$less[sel.l, "q.val"])
  if (nrow(df_downregulated) > 1){
    df_downregulated$number <- 1 : nrow(df_downregulated)
  }
  write.csv(df_downregulated, paste0("./gage_pathview/", "gage_", sample_list, "_upreg.csv"))
  pv_downregulated <- sapply(path.ids.l, function(pid) pathview(gene.data = deseq2.fc,
                                                                pathway.id = pid,
                                                                species = "mmu",
                                                                limit = list(gene = 3, cpd = 2),
                                                                bins = list(gene = 20, cpd = 20),
                                                                out.suffix = paste0(sample_list, "_downreg"),
                                                                kegg.dir = "./gage_pathview/"))
  
  
})







