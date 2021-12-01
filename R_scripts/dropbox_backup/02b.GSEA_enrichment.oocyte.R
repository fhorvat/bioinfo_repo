### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/review/GO_terms/GSEA_enrichment")

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
experiment <- "hamster_oocyte_Mov10l.RNAseq"
experiment_name <- "hamster_oocyte_Mov10l.RNAseq"

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
analysis_path <- file.path(base_path, "Analysis/expression.added_PIWIL3")


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
diffExp_path <- list.files(diffExp_path, "protein_coding.diffExp.DESeq2.genotype.significant_results.xlsx", full.names = T)


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
comparison_list <- c("Mov10l_KO_vs_Mov10l_WT", "Mov10l_HET_vs_Mov10l_WT")

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
  # do GO enrichment across 3 GO categories
  go_categories <- c("CC", "MF", "BP")
  
  # enrich
  go_enrich <- purrr::map(go_categories, function(go_cat){
    
    # enrich GO terms
    ego <- clusterProfiler::gseGO(gene = geneList,
                                  keyType = "ENSEMBL",
                                  OrgDb = org.Mm.eg.db,
                                  ont = go_cat,
                                  minGSSize = 3,
                                  maxGSSize = 800,
                                  pvalueCutoff = 0.05,
                                  verbose = TRUE)
    
    if(nrow(ego) > 0){
      
      # save to pdf
      pdf(file = file.path(outpath, str_c("GSEAEnrich.visualization", comparison, go_cat, "full.pdf", sep = ".")),
          width = 15,
          height = 10)
      tryCatch({print(dotplot(ego, showCategory = 10, split = ".sign") + facet_grid( . ~ .sign))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      tryCatch({print(emapplot(ego, showCategory = 10))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      tryCatch({print(cnetplot(ego, categorySize = "pvalue", foldChange = geneList, showCategory = 3))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      tryCatch({print(ridgeplot(ego) + labs(x = "enrichment distribution"))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      dev.off()
      
      # simplify enriched GO terms by removing redundancy
      ego_simple <- clusterProfiler::simplify(x = ego,
                                              cutoff = 0.7,
                                              by = "p.adjust",
                                              select_fun = min)
      
      # save visualizations
      if(nrow(ego_simple) > 0){
        
        # save to pdf
        pdf(file = file.path(outpath, str_c("GSEAEnrich.visualization", comparison, go_cat, "simple.pdf", sep = ".")),
            width = 15,
            height = 10)
        tryCatch({print(dotplot(ego_simple, showCategory = 10, split = ".sign") + facet_grid( . ~ .sign))}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
        tryCatch({print(emapplot(ego_simple, showCategory = 10))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
        tryCatch({print(cnetplot(ego_simple, categorySize = "pvalue", foldChange = geneList, showCategory = 3))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
        tryCatch({print(ridgeplot(ego_simple) + labs(x = "enrichment distribution"))}, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
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
  set_names(., comparison_list)

# save results of enrichment
saveRDS(go_enrich_list, file = file.path(outpath, str_c(experiment_name, 
                                                        "DESeq2.protein_coding.significant_results.GSEA.RDS", sep = ".")))


### write results
if(length(unlist(go_enrich_list)) > 0){
  
  # loop through comparison
  purrr::map(comparison_list, function(comparison){
    
    # get enrichment table for one stage
    go_enrich <- go_enrich_list[[comparison]]
    
    # join results - full GO enrich
    ego_full_df <- 
      purrr::map(go_enrich, function(x) x$full_ego) %>% 
      dplyr::bind_rows(.)
    
    # continue if not empty
    if(nrow(ego_full_df) > 0){
      
      # clean and save table
      ego_full_df %>% 
        dplyr::filter(!duplicated(ID)) %>% 
        dplyr::rename(GO_id = ID) %>% 
        dplyr::select(category = GO_category, GO_id, Term = Description, adj_pval = p.adjust, Genes = core_enrichment) %>% 
        dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                      Genes = str_replace_all(Genes, "/", ", ")) %T>% 
        readr::write_csv(., file = file.path(outpath, str_c(experiment, "GSEAEnrich.results", comparison, "full.csv", sep = ".")))
      
      # simplified GO enrich
      ego_simple_df <- 
        purrr::map(go_enrich, function(x) x$simple_ego) %>% 
        dplyr::bind_rows(.) %>% 
        dplyr::filter(!duplicated(ID)) %>% 
        dplyr::rename(GO_id = ID) %>%
        dplyr::select(category = GO_category, GO_id, Term = Description, adj_pval = p.adjust, Genes = core_enrichment) %>% 
        dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                      Genes = str_replace_all(Genes, "/", ", ")) %T>% 
        readr::write_csv(., file = file.path(outpath, str_c(experiment, "GSEAEnrich.results", comparison, "simple.csv", sep = ".")))
      
    }
    
  })
  
}


