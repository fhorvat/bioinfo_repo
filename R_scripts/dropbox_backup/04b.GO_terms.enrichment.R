### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()
# LD_LIBRARY_PATH=/common/WORK/fhorvat/programi/udunits2/udunits-2.2.26/lib:$LD_LIBRARY_PATH R
# install.packages(pkgs = "udunits2",
#                  configure.args = "--with-udunits2-lib=/common/WORK/fhorvat/programi/udunits2/udunits-2.2.26/lib --with-udunits2-include=/common/WORK/fhorvat/programi/udunits2/udunits-2.2.26/include")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_expression/GO_terms")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# modified GOplot::GOBubble which uses ggrepel::geom_text_repel instead of ggplot2::geom_text
GOBubble2 <- function(data, display, title, colour, labels, ID = T, table.legend = T, table.col = T, bg.col = F, pvalue_threshold = 0.01){
  
  # ID = T
  # table.legend = T
  # table.col = T
  # bg.col = F
  # pvalue_threshold = 0.01
  # data = circ
  # labels = 2
  # display = "multiple"
  # bg.col = T
  # title = str_c("Enriched GO terms", regulation, "genes", sep = " ")
  
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

# path to differentialy expressed genes
diff_path <- file.path(inpath, "DE_genes.GO_terms.csv")

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# GO terms info path
go_info_path <- file.path(genome_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.GOterms.csv")

# gene info path
gene_info_path <- file.path(genome_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv")

######################################################## READ DATA
# read list of diff. expressed genes
diff_df <- readr::read_csv(diff_path)

# read gene info
go_info <- readr::read_csv(go_info_path)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

######################################################## MAIN CODE
# prepare gene list
## feature 1: numeric vector
geneList <- diff_df$log2FoldChange_NW
names(geneList) <- diff_df$gene_id
geneList <- sort(geneList, decreasing = TRUE)

### GO terms enrichment
# do GO enrichment across 3 GO categories
go_categories <- c("CC", "MF", "BP")

go_enrich <- purrr::map(go_categories[1], function(go_cat){

  # enrich GO terms
  ego <- clusterProfiler::enrichGO(gene = diff_df$gene_id,
                                   universe = gene_info$gene_id,
                                   keyType = "ENSEMBL",
                                   OrgDb = org.Mm.eg.db,
                                   ont = go_cat,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)

  if(nrow(ego) > 0){

    # # save to pdf
    # pdf(file = file.path(outpath, str_c("GOEnrich.visualization", go_cat, "full.pdf", sep = ".")),
    #     width = 15,
    #     height = 10)
    # tryCatch({print(barplot(ego, drop = TRUE, showCategory = 15))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    # tryCatch({print(dotplot(ego))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    # tryCatch({print(emapplot(ego))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    # tryCatch({print(goplot(ego, geom = "label"))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    # tryCatch({print(cnetplot(ego, categorySize = "pvalue", foldChange = geneList))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    # dev.off()

    # simplify enriched GO terms by removing redundancy
    ego_simple <- clusterProfiler::simplify(x = ego,
                                            cutoff = 0.7,
                                            by = "p.adjust",
                                            select_fun = min)

    # save visualizations
    if(nrow(ego_simple) > 0){

      # # save to pdf
      # pdf(file = file.path(outpath, str_c("GOEnrich.visualization", go_cat, "simple.pdf", sep = ".")),
      #     width = 15,
      #     height = 10)
      # tryCatch({print(barplot(ego_simple, drop = TRUE, showCategory = 15))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      # tryCatch({print(dotplot(ego_simple))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      # tryCatch({print(emapplot(ego_simple))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      # tryCatch({print(goplot(ego_simple, geom = "label"))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      # tryCatch({print(cnetplot(ego_simple, categorySize = "pvalue", foldChange = geneList))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      # dev.off()

      # convert to tibble
      ego_df_simple <-
        ego_simple %>%
        as.data.frame(.) %>%
        as.tibble(.) %>%
        dplyr::mutate(GO_category = go_cat) %>%
        dplyr::select(ID, GO_category, everything())

    }

    # convert to tibble
    ego_df <-
      ego %>%
      as.data.frame(.) %>%
      as.tibble(.) %>%
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

# save results of enrichment
saveRDS(go_enrich, file = file.path(outpath, "go_enrich.RDS"))

# read enrichment calculations
go_enrich <- readRDS(file = file.path(outpath, "go_enrich.RDS"))

### write results, plot bubble plot
# join results - full GO enrich
ego_full_df <- 
  purrr::map(go_enrich, function(x) x$full_ego) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::filter(!duplicated(ID)) %>% 
  dplyr::rename(GO_id = ID) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("GOEnrich.results.full.csv"))) %>% 
  dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
  dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                Genes = str_replace_all(Genes, "/", ", "))

# simplified GO enrich
ego_simple_df <- 
  purrr::map(go_enrich, function(x) x$simple_ego) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::filter(!duplicated(ID)) %>% 
  dplyr::rename(GO_id = ID) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("GOEnrich.results.simplified.csv"))) %>% 
  dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
  dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                Genes = str_replace_all(Genes, "/", ", "))

# rename columns in results
genes_info_df <-
  diff_df %>% 
  dplyr::rename(padj = padj_NW, ID = gene_name, logFC = log2FoldChange_NW)

### plot results
for(result in c("full", "simple")){
  
  # set results table
  if(result == "full"){
    circ <- GOplot::circle_dat(terms = ego_full_df, genes = genes_info_df)
  }else{
    circ <- GOplot::circle_dat(terms = ego_simple_df, genes = genes_info_df)
  }
  
  # GO bubble plot with table
  png(file = file.path(outpath, str_c("GOEnrich.bubble.single.table.", result, ".png")),
      width = 1800, 
      height = 1000, 
      units = "px", 
      type = "cairo")
  GOBubble2(data = circ,
            labels = 1,
            display = "single",
            table.legend = T,
            pvalue_threshold = 0.01,
            title = str_c("Enriched GO terms - ", result))
  dev.off()
  
  # GO bubble plot facet
  png(file = file.path(outpath, str_c("GOEnrich.bubble.facet.", result, ".png")),
      width = 1000, 
      height = 1000, 
      units = "px", 
      type = "cairo")
  GOBubble2(data = circ,
            labels = 1,
            display = "multiple",
            bg.col = T,
            pvalue_threshold = 0.01,
            title = str_c("Enriched GO terms - ", result))
  dev.off()
  
}

