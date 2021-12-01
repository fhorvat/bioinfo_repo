### INFO: CNOT6L GO terms enrichment analysis
### DATE: Thu Apr 05 22:09:25 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(clusterProfiler)
library(GOplot)
library(org.Mm.eg.db)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.geneInfo.csv"

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_coding <- 
  ensembl_genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# get gene symbols from org.Mm.eg.db
protein_coding_annotation <- clusterProfiler::bitr(protein_coding, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Mm.eg.db")

### GO terms enrichment - enrich and write data
# loop through stages
for(filt_stage in c("GV", "MII", "1C")){
  
  # read results 
  results_df <- 
    read_csv(file = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv"))) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1)) %>% 
    left_join(., protein_coding_annotation, by = c("gene_id" = "ENSEMBL")) %>% 
    dplyr::select(gene_id, gene_name = SYMBOL, entrezID = ENTREZID, log2FoldChange, padj) %>% 
    dplyr::filter(!is.na(entrezID))
  
  ## separate GO enrichment for up- and down-regulated genes
  for(regulation in c("upregulated", "downregulated")){
    
    # filter by regulation
    if(regulation == "upregulated"){
      results_filt <- dplyr::filter(results_df, log2FoldChange > 0, padj < 0.1)
    }else{
      results_filt <- dplyr::filter(results_df, log2FoldChange < 0, padj < 0.1)
    }
    
    # GO enrichment and visualization
    if(nrow(results_filt) > 0){
      
      ### GO enrichment with clusterProfile
      # do GO enrichment across 3 GO categories
      go_enrich <- 
        lapply(c("CC", "MF", "BP"), function(go_cat){
          
          # enrich GO terms
          ego <- 
            clusterProfiler::enrichGO(gene = results_filt$entrezID,
                                      universe = protein_coding_annotation$ENTREZID,
                                      OrgDb = org.Mm.eg.db,
                                      ont = go_cat,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.01,
                                      qvalueCutoff = 0.05,
                                      readable = TRUE)
          
          if(nrow(ego) > 0){
            
            # saveGO enriched map
            png(file = file.path(outpath, "results", str_c("GOEnrich.map.CNOT6L.", filt_stage, ".", regulation, ".", go_cat, ".whole.png")),
                width = 1000, height = 1000, units = "px", type = "cairo")
            enrichMap(ego)
            dev.off()
            
          }
          
          # simplify enriched GO terms by removing redundancy 
          ego_simpl <- 
            clusterProfiler::simplify(x = ego,
                                      cutoff = 0.7,
                                      by = "p.adjust",
                                      select_fun = min)
          
          # save simplified GO enriched map
          if(nrow(ego_simpl) > 0){
            
            png(file = file.path(outpath, "results", str_c("GOEnrich.map.CNOT6L.", filt_stage, ".", regulation, ".", go_cat, ".simplified.png")),
                width = 1000, height = 1000, units = "px", type = "cairo")
            enrichMap(ego_simpl)
            dev.off()
            
          }
          
          # convert to data.frame
          ego_df <- 
            ego_simpl %>% 
            as.data.frame(.) %>% 
            as.tibble(.) %>% 
            dplyr::mutate(GO_category = go_cat) %>% 
            dplyr::select(ID, GO_category, everything())
          
          # return
          return(ego_df)
          
        }) %>% 
        dplyr::bind_rows(.)
      
      
      ### save and visualize results
      if(nrow(go_enrich) > 0){
        
        # write results
        go_enrich %<>% 
          dplyr::mutate(geneID = str_replace_all(geneID, "/", ", ")) %>% 
          dplyr::filter(!duplicated(ID)) %>% 
          dplyr::rename(GO_id = ID) %T>%
          readr::write_csv(., path = file.path(outpath, "results", str_c("GOEnrich.results.CNOT6L.", filt_stage, ".", regulation, ".csv")))
        
      }
    }
  }
}

### visualize GO terms enrichment using GOplot
# loop through stages
for(filt_stage in c("GV", "MII", "1C")){
  
  # loop through regulation
  for(regulation in c("upregulated", "downregulated")){
    
    # check if file exists
    if(file.exists(file.path(outpath, "results", str_c("GOEnrich.results.CNOT6L.", filt_stage, ".", regulation, ".csv")))){
      
      ### prepare data - GO enrichment data.frame
      # read data
      go_enrich_df <- 
        readr::read_csv(file = file.path(outpath, "results", str_c("GOEnrich.results.CNOT6L.", filt_stage, ".", regulation, ".csv"))) %>% 
        dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
        dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")))
      
      ### prepare data - genes info data.frame
      # read results 
      genes_info_df <- 
        read_csv(file = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv"))) %>% 
        dplyr::mutate(padj = replace(padj, is.na(padj), 1)) %>% 
        left_join(., protein_coding_annotation, by = c("gene_id" = "ENSEMBL")) %>% 
        dplyr::filter(!is.na(ENTREZID)) %>% 
        dplyr::rename(ID = SYMBOL, logFC = log2FoldChange)
      
      # filter by regulation
      if(regulation == "upregulated"){
        genes_info_df <- dplyr::filter(genes_info_df, logFC > 0, padj < 0.1)
      }else{
        genes_info_df <- dplyr::filter(genes_info_df, logFC < 0, padj < 0.1)
      }
      
      # generate the plotting object
      circ <- GOplot::circle_dat(terms = go_enrich_df, genes = genes_info_df)
      
      # GO bubble plot without table
      png(file = file.path(outpath, "results", str_c("GOEnrich.bubble.single.CNOT6L.", filt_stage, ".", regulation, ".png")),
          width = 1000, height = 1000, units = "px", type = "cairo")
      GOBubble2(data = circ,
                labels = 2,
                display = "single",
                table.legend = F,
                title = str_c("Enriched GO terms - CNOT6L", filt_stage, regulation, "genes", sep = " "))
      dev.off()
      
      # GO bubble plot with table
      png(file = file.path(outpath, "results", str_c("GOEnrich.bubble.single.table.CNOT6L.", filt_stage, ".", regulation, ".png")),
          width = 1800, height = 1000, units = "px", type = "cairo")
      GOBubble2(data = circ,
                labels = 2,
                display = "single",
                table.legend = T,
                title = str_c("Enriched GO terms - CNOT6L", filt_stage, regulation, "genes", sep = " "))
      dev.off()
      
      # GO bubble plot facet
      png(file = file.path(outpath, "results", str_c("GOEnrich.bubble.facet.CNOT6L.", filt_stage, ".", regulation, ".png")),
          width = 1000, height = 1000, units = "px", type = "cairo")
      GOBubble2(data = circ,
                labels = 2,
                display = "multiple",
                bg.col = T,
                title = str_c("Enriched GO terms - CNOT6L", filt_stage, regulation, "genes", sep = " "))
      dev.off()
      
    }
  }
}
