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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Mm.eg.db)
library(openxlsx)
library(gprofiler2)

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
  
  # do the enrichment
  gostres <- gost(query = names(geneList), 
                  organism = "mmusculus", ordered_query = TRUE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL, 
                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
}) %>% 
  set_names(., comparison_list)

# save results of enrichment
saveRDS(go_enrich_list, file = file.path(outpath, str_c(experiment_name, 
                                                        "DESeq2.protein_coding.significant_results.GO_terms.RDS", sep = ".")))


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
        dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
        dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                      Genes = str_replace_all(Genes, "/", ", ")) %T>% 
        readr::write_csv(., file = file.path(outpath, str_c(experiment, "GOEnrich.results", comparison, "full.csv", sep = ".")))
      
      # simplified GO enrich
      ego_simple_df <- 
        purrr::map(go_enrich, function(x) x$simple_ego) %>% 
        dplyr::bind_rows(.) %>% 
        dplyr::filter(!duplicated(ID)) %>% 
        dplyr::rename(GO_id = ID) %>%
        dplyr::select(category = GO_category, ID = GO_id, Term = Description, Genes = geneID, adj_pval = p.adjust) %>% 
        dplyr::mutate(category = factor(category, levels = c("CC", "MF", "BP")), 
                      Genes = str_replace_all(Genes, "/", ", ")) %T>% 
        readr::write_csv(., file = file.path(outpath, str_c(experiment, "GOEnrich.results", comparison, "simple.csv", sep = ".")))
      
    }
    
  })
  
}




### pathview
# additional libraries
library(biomaRt)
library(pathview)
library(gage)
data(kegg.gs)

# KEGG set for mouse
kg.mouse <- kegg.gsets("mouse")
kegg.gs <- kg.mouse$kg.sets[kg.mouse$sigmet.idx]

# get info about genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", 
                host = "http://jan2020.archive.ensembl.org")

# get entrez IDs
entrez_info <-
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart) %>%
  as_tibble(.) %>% 
  dplyr::rename(mmusculus_homolog_ensembl_gene = ensembl_gene_id,
                mmusculus_homolog_entrez_id = entrezgene_id)

### write results, plot 
purrr::map(comparison_list, function(comparison){
  
  # get results table for one comparison, add mouse homologs
  diff_df <- 
    diffExp_results[[comparison]] %>% 
    dplyr::left_join(., mouse_homologs, by = "gene_id") %>% 
    dplyr::filter(!is.na(mmusculus_homolog_ensembl_gene)) %>% 
    dplyr::left_join(., entrez_info, by = "mmusculus_homolog_ensembl_gene") %>% 
    dplyr::filter(!is.na(mmusculus_homolog_entrez_id))
  
  # prepare gene list
  geneList <- diff_df$log2FoldChange
  names(geneList) <- diff_df$mmusculus_homolog_entrez_id
  geneList <- sort(geneList, decreasing = TRUE)
  
  # gage
  fc.kegg.p <- gage(geneList, gsets = kegg.gs, ref = NULL, samp = NULL)
  
  ### get upregulated pathways
  # upregulated pathways
  sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
  path.ids <- rownames(fc.kegg.p$greater)[sel]
  path.ids <- str_remove(path.ids, " .*")
   
  ## save and plot
  if(length(path.ids) > 1){
    
    # save as table
    df_upregulated <- 
      tibble(id = rownames(fc.kegg.p$greater)[sel],
             q.value = fc.kegg.p$greater[sel, "q.val"]) %>% 
      dplyr::mutate(number = 1:nrow(.)) %T>% 
      readr:::write_csv(., file.path(outpath, "gage_pathview", str_c("gage.", comparison, ".upregulated.csv")))
    
    # plot
    pv_upregulated <- sapply(path.ids, function(pid) pathview(gene.data = geneList,
                                                              pathway.id = pid,
                                                              species = "mmu",
                                                              limit = list(gene = 3, cpd = 2),
                                                              bins = list(gene = 20, cpd = 20),
                                                              out.suffix = str_c(comparison, "_upregulated"),
                                                              kegg.dir = file.path(outpath, "gage_pathview")))
    
  }

  
  ### get downregulated pathways
  # downregulated pathways
  sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
  path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
  path.ids.l <- str_remove(path.ids.l, " .*")
  
  ## save and plot
  if(length(path.ids.l) > 1){
    
    # save as table
    df_upregulated <- 
      tibble(id = rownames(fc.kegg.p$less)[sel.l],
             q.value = fc.kegg.p$less[sel.l, "q.val"]) %>% 
      dplyr::mutate(number = 1:nrow(.)) %T>% 
      readr:::write_csv(., file.path(outpath, "gage_pathview", str_c("gage.", comparison, ".downregulated.csv")))
    
    # plot
    pv_upregulated <- sapply(path.ids.l, function(pid) pathview(gene.data = geneList,
                                                              pathway.id = pid,
                                                              species = "mmu",
                                                              limit = list(gene = 3, cpd = 2),
                                                              bins = list(gene = 20, cpd = 20),
                                                              out.suffix = str_c(comparison, "_downregulated"),
                                                              kegg.dir = file.path(outpath, "gage_pathview")))
    
  }
  
})







