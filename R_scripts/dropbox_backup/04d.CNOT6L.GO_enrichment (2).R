### INFO: CNOT6L GO terms enrichment analysis
### DATE: Thu Apr 05 22:09:25 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

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
library(biomaRt)
library(goseq)
library(GO.db)
library(clusterProfiler)
library(biomaRt)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_coding <- 
  ensembl_genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

### GO terms enrichment
# loop through stages
for(filt_stage in unique(sample_table_dds$stage)){
  
  # read results 
  results_df <- 
    read_csv(file = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv"))) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1))
  
  ### clusterProfile
  # get upregulated genes 
  genes <- 
    results_df %>% 
    dplyr::filter(padj < 0.1, log2FoldChange > 0) %$% 
    gene_id
  
  # do GO enrichment across 3 GO categories
  go_enrich <- lapply(c("CC", "MF", "BP"), function(go_cat){
    
    # enrich and simplify GO terms by removing redundancy 
    ego <- 
      clusterProfiler::enrichGO(gene = genes,
                                universe = protein_coding,
                                keyType = "ENSEMBL", 
                                OrgDb = org.Mm.eg.db,
                                ont = go_cat,
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2,
                                readable = TRUE) %>% 
      clusterProfiler::simplify(x = ., 
                                cutoff = 0.7, 
                                by = "p.adjust", 
                                select_fun = min)
    
    # save simplified GO enriched map
    png(file = file.path(outpath, "results", str_c("GOEnrichMap.", go_cat, ".CNOT6L.", filt_stage, ".upregulated.png")),
        width = 1000, height = 1000, units = "px", type = "cairo")
    enrichMap(ego)
    dev.off()
    
    # convert to data.frame
    ego_df <- 
      ego %>% 
      as.data.frame(.) %>% 
      as.tibble(.) %>% 
      dplyr::mutate(GO_category = go_cat) %>% 
      dplyr::select(ID, GO_category, everything())
    
    # return
    return(ego_df)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  ## Venn diagram
  png(file = file.path(outpath, "results", "Venn.Ma_vs_CNOT6L.significant.upregulated.png"),
      width = 1000, height = 1000, units = "px", type = "cairo")
  up_plot <- venn.diagram(x = list(Ma = na.omit(results_Ma$gene_id),
                                   CNOT6L = results %>% dplyr::filter(stage == signif_stage, signif == "yes") %$% gene_id),
                          filename = NULL,
                          fill = c("blue", "red"),
                          alpha = c(0.5, 0.5),
                          cex = 4,
                          cat.cex = 2,
                          main = "MII KO vs. WT upregulated genes",
                          main.cex = 2)
  grid.draw(up_plot)
  dev.off()
  
  # ### goseq
  # # create named genes vector
  # genes <- 
  #   as.integer(results_df$padj < 0.1 & results_df$log2FoldChange > 0) %>% 
  #   magrittr::set_names(., results_df$gene_id)
  # 
  # # fit the Probability Weighting Function (PWF)
  # pwf <- goseq::nullp(genes, "mm10", "ensGene")
  # 
  # # Wallenius approximation (GO:CC, GO:BP, GO:MF)
  # GO.wall <- goseq::goseq(pwf, "mm10", "ensGene", test.cats = "GO:CC")
  # 
  # # get enriched GO terms
  # GO.wall.enriched <- GO.wall[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05, ]
  
  # ### GAGE - BP, CC, MF
  # ## get additional info about genes from Biomart
  # # load Mart of mouse database from ensembl
  # # sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
  # mart <- "error"
  # count <- 0
  # while(class(mart) == "character"){
  #   count <- count + 1
  #   print(str_c("mmusculus_gene_ensembl", " ", count))
  #   mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"), 
  #                    error = function(x) return("error"))
  # }
  # 
  # # get entrezgene ID
  # ensembl_info <- 
  #   getBM(attributes = c("ensembl_gene_id", "entrezgene"), mart = mart) %>% 
  #   as.tibble(.) %>% 
  #   dplyr::rename(gene_id = ensembl_gene_id)
  # 
  # # create gage gset for mouse
  # go_mouse <- gage::go.gsets(species = "mouse")
  # 
  # # join with entrez ID
  # results_df_GO <- 
  #   results_df %>% 
  #   dplyr::left_join(., ensembl_info, by = "gene_id") 
  # 
  # # create vector of named log2FC changes (names = entrez IDs)
  # logFC_vector <-
  #   results_df_GO$log2FoldChange %>%
  #   set_names(., results_df_GO$entrezgene)
  # 
  # # run GAGE
  # go_gage <- gage::gage(exprs = logFC_vector, gsets = go_mouse$go.sets[go_mouse$go.subs$MF], ref = NULL, samp = NULL)
  # sel <- go_gage$greater[, "q.val"] < 0.1 & !is.na(go_gage$greater[, "q.val"])
  # path.ids <- rownames(go_gage$greater)[sel]
  # sel.l <- go_gage$less[, "q.val"] < 0.1 & !is.na(go_gage$less[,"q.val"])
  # path.ids.l <- rownames(go_gage$less)[sel.l]
  
  
}
