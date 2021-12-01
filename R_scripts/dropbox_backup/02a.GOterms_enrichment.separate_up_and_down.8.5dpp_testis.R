### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/review/GO_terms")

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

######################################################## PATH VARIABLES
### experiment
# set experiment name
experiment <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"
experiment_name <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"

### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### other experiment paths
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", experiment)

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mesAur1")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded")


### documentation
# set ensembl version
ensembl_version <- 99

# sample table path
sample_table_path <- list.files(documentation_path, ".*\\.sampleTable\\.csv", full.names = T)

# stats and tracks path
stats_and_tracks_path <- list.files(path = mapped_path, pattern = ".*\\.stats_and_tracks\\.csv", full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.FPKM_mean\\.csv$"), full.names = T)

# differential expression results path
diffExp_path <- file.path(analysis_path, "results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq")
diffExp_path <- list.files(diffExp_path, "protein_coding.diffExp.DESeq2.genotype_age.significant_results.xlsx", full.names = T)


### genome
# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# gene info path
genes_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# GO terms info
go_terms_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.GOterms.csv$"), full.names = T)

# mouse homologs path
mouse_homologs_path <- list.files(genome_path, pattern = str_c("ensembl.", ensembl_version, ".*mmusculus_homologs.csv$"), full.names = T)

# mouse genes info
ensembl_version <- 99
genes_info_path_mouse <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"
genes_info_path_mouse <- list.files(path = genes_info_path_mouse, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

######################################################## READ DATA
# read gene info
genes_info <- 
  readr::read_csv(genes_info_path) %>% 
  dplyr::filter(gene_biotype == "protein_coding")

# read GO terms info
go_term_info <- readr::read_csv(go_terms_info_path)

# read mouse homologs
mouse_homologs <- readr::read_csv(mouse_homologs_path)

# read mouse genes info
genes_info_mouse <- 
  readr::read_csv(genes_info_path_mouse) %>% 
  dplyr::filter(gene_biotype == "protein_coding")


### read list of diff. expressed genes
# set comparison
comparison_list <- c("Mov10l1_KO_8.5dpp_vs_Mov10l1_WT")

# read different comparisons
diffExp_results <- purrr::map(comparison_list, function(comparison){
  
  # read sheet
  read.xlsx(xlsxFile = diffExp_path, sheet = comparison) %>% 
    as_tibble(.)
  
}) %>% 
  set_names(comparison_list)

######################################################## MAIN CODE
### do the enrichment for all comparisons
go_enrich_list <- purrr::map(comparison_list, function(comparison){
  
  # get results table for one stage, add mouse homologs
  diff_df <- 
    diffExp_results[[comparison]] %>% 
    left_join(., mouse_homologs, by = "gene_id") %>% 
    dplyr::filter(!is.na(mmusculus_homolog_ensembl_gene)) %>% 
    dplyr::filter(!(mmusculus_homolog_ensembl_gene %in% mmusculus_homolog_ensembl_gene[duplicated(mmusculus_homolog_ensembl_gene)]))
  
  # prepare gene list
  geneList <- diff_df$log2FoldChange
  names(geneList) <- diff_df$mmusculus_homolog_ensembl_gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  ### GO terms enrichment
  # do enrichment separately for upregulated and downregulated genes
  go_enrich_regulation <- purrr::map(c("upregulated", "downregulated"), function(regulation){
    
    # filter based on regulation
    if(regulation == "upregulated"){
      geneList <- geneList[geneList > 0]
    }else{
      if(regulation == "downregulated"){
        geneList <- geneList[geneList < 0]
      }
    }
    
    # do GO enrichment across 3 GO categories
    go_categories <- c("CC", "MF", "BP")
    
    # enrich
    go_enrich <- purrr::map(go_categories, function(go_cat){
      
      # enrich GO terms
      ego <- clusterProfiler::enrichGO(gene = names(geneList),
                                       universe = unique(na.omit(mouse_homologs$mmusculus_homolog_ensembl_gene)),
                                       keyType = "ENSEMBL",
                                       OrgDb = org.Mm.eg.db,
                                       ont = go_cat,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.01,
                                       qvalueCutoff = 0.05,
                                       readable = TRUE)
      
      if(nrow(ego) > 0){
        
        # save to pdf
        pdf(file = file.path(outpath, str_c("GOEnrich.visualization", comparison, go_cat, regulation, "full.pdf", sep = ".")),
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
          pdf(file = file.path(outpath, str_c("GOEnrich.visualization", comparison, go_cat, regulation, "simple.pdf", sep = ".")),
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
      
    }) %>% 
      set_names(go_categories)
    
  }) %>% 
    set_names(c("upregulated", "downregulated"))
  
}) %>% 
  set_names(., comparison_list)

# save results of enrichment
saveRDS(go_enrich_list, file = file.path(outpath, str_c(experiment_name, 
                                                        "separate_regulation", 
                                                        "DESeq2.protein_coding.significant_results.GO_terms.RDS", sep = ".")))


### write results
if(length(unlist(go_enrich_list)) > 0){
  
  # loop through comparison
  purrr::map(comparison_list, function(comparison){
    
    # get enrichment table for one stage
    go_enrich <- go_enrich_list[[comparison]]
    
    # loop through regulation
    purrr::map(c("upregulated", "downregulated"), function(regulation){
      
      # join results - full GO enrich
      ego_full_df <- 
        go_enrich[[regulation]] %>% 
        purrr::map(., function(x) x$full_ego) %>% 
        dplyr::bind_rows(.)
      
      # continue if not empty
      if(nrow(ego_full_df) > 0){
        
        # clean and save table
        ego_full_df %>% 
          dplyr::filter(!duplicated(ID)) %>% 
          dplyr::rename(GO_id = ID) %>% 
          dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
          dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                        Genes = str_replace_all(Genes, "/", ", ")) %T>% 
          readr::write_csv(., file = file.path(outpath, str_c(experiment, "GOEnrich.results", comparison, regulation, "full.csv", sep = ".")))
        
        # simplified GO enrich
        ego_simple_df <- 
          purrr::map(go_enrich[[regulation]], function(x) x$simple_ego) %>% 
          dplyr::bind_rows(.) %>% 
          dplyr::filter(!duplicated(ID)) %>% 
          dplyr::rename(GO_id = ID) %>%
          dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
          dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                        Genes = str_replace_all(Genes, "/", ", ")) %T>% 
          readr::write_csv(., file = file.path(outpath, str_c(experiment, "GOEnrich.results", comparison, regulation, "simple.csv", sep = ".")))
        
      }
      
      # return 
      return(regulation)
      
    })
    
    # return
    return(comparison)
    
  })
  
}


