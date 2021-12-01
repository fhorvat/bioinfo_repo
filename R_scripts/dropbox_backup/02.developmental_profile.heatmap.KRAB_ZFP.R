### INFO: 
### DATE: Tue Jun 04 09:35:24 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/KRAB_ZFP")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# changes draw_legend function in pheatmap
my.draw.legend <-  function (color, breaks, legend, ...){
  
  library(grid)
  
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos = height * legend_pos + (unit(1, "npc") - height)
  breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
  breaks = height * breaks + (unit(1, "npc") - height)
  h = breaks[-1] - breaks[-length(breaks)]
  h2 = breaks[length(breaks)] - breaks[1]
  rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = NA))
  rect1 = rectGrob(x = 0, y = breaks[1], width = unit(10, "bigpts"), height = h2, hjust = 0, vjust = 0, gp = gpar(fill = NA, col = "grey60"))
  text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
  res = grobTree(rect, text, rect1)
  
  return(res)
  
}

# change pheatmap:::draw_legendwith my.draw.legend
assignInNamespace("draw_legend", my.draw.legend, ns = "pheatmap")

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)


### experiments
# experiment list
experiment_list <- c("Deng_2014_Science_GSE45719", 
                     "Veselovska_2015_GenomeBiol_GSE70116",
                     "Yamaguchi_2013_CellRes_GSE41908", 
                     "Gan_2013_NatCommun_GSE35005", 
                     "ENCODE_2014_Nature_GSE49417")

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/developmental_profile_expression")

# documentation path
documentation_path <- file.path(base_path, "Documentation")

# summarizedExperiment path
se_main_path <- file.path(base_path, "summarizedExperiments")

# FPKM path
fpkm_main_path <- file.path(base_path, "FPKMs")


### documentation and data
# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable.csv", full.names = T)

# fpkm paths
fpkm_paths <- 
  list.files(path = fpkm_main_path, pattern = str_c("ensembl.", ensembl_version, ".*\\.FPKM_statistics\\.csv"), full.names = T) %>% 
  .[str_detect(., str_c(experiment_list, collapse = "|"))]

# KRAB genes path
krab_genes_path <- file.path(inpath, "ensembl.93.KRAB_domain.pfam01352.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read sample table
sample_table <- 
  data.table::fread(sample_table_path) %>% 
  .[experiment %in% experiment_list]

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read FPKMs
fpkm_tb <- 
  purrr::map(fpkm_paths, function(path){
    
    # read 
    data.table::fread(path) %>% 
      .[, experiment := str_extract(path, str_c(experiment_list, collapse = "|"))] %>% 
      .[]
    
  }) %>% 
  rbindlist(.) 

# read list of genes with KRAB domain
krab_genes <- read_csv(krab_genes_path)

######################################################## MAIN CODE
# clean stage names
stage_clean <- 
  tibble(stage = c("PGC_9.5", "PGC_11.5", "PGC_13.5_m", "PGC_13.5_f", 
                   "primitive_SG_A", "SG_B", "leptotene_SC", "pachytene_SC", "round_ST", 
                   "nonGrowing_oocytes", "growing_oocytes_d8_14", "GV_oocytes",
                   "zygote", "2C", "4C", "8C", "16C", "mid_blast",
                   "placenta_adult8wks"), 
         stage_clean = c("PGC 9.5", "PGC 11.5", "male PGC 13.5", "female PGC 13.5", 
                         "primitive SG A", "SG B", "leptotene SC", "pachytene SC", "round ST",
                         "non-growing oocyte", "growing oocyte days 8-14", "GV ",
                         "zygote", "2-cell", "4-cell", "8-cell", "16-cell", "blastocyst 94h", 
                         "placenta")) %>%
  mutate(stage_clean = factor(stage_clean, levels = stage_clean))

# filter FPKM table to include only gene of interest
fpkm_tb_filtered <-   
  fpkm_tb[(gene_id %in% krab_genes$gene_id) & (stage %in% stage_clean$stage), ] %>% 
  as_tibble(.) %>% 
  left_join(., stage_clean, by = "stage")

# save as wide table
invisible(fpkm_tb_filtered %>% 
            dplyr::select(gene_id, stage_clean, avg_fpkm) %>%
            tidyr::spread(key = stage_clean, value = avg_fpkm) %>% 
            dplyr::left_join(., genes_info, by = "gene_id") %>% 
            dplyr::arrange(desc(`non-growing oocyte`)) %>% 
            dplyr::select(gene_id, gene_name, everything()) %T>%
            write_csv(., file.path(outpath, str_c("FPKM_average.developmental_profile.ensembl", ensembl_version, "KRAB_domain.pfam01352.csv", sep = "."))))

# get selected genes' FPKMs, create matrix for heatmap
fpkm_heatmap_tb <- 
  fpkm_tb_filtered %>%
  mutate(log10_fpkm = log10(avg_fpkm + 0.001)) %>% 
  dplyr::left_join(., genes_info, by = "gene_id") %>% 
  dplyr::mutate(gene_id = str_c(gene_name, " (", gene_id, ")")) %>% 
  dplyr::select(gene_id, stage_clean, log10_fpkm) %>%
  tidyr::spread(key = stage_clean, value = log10_fpkm) %>% 
  dplyr::arrange(desc(`non-growing oocyte`)) %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_id") %>% 
  as.matrix(.)

### plot heatmap
# get columns and rows to put gaps
which.gap.col <- which(colnames(fpkm_heatmap_tb) %in% c("female PGC 13.5", "round ST", "fully-grown oocyte", "GV ", "blastocyst 94h", "placenta"))

# plot heatmap with annotation
pheatmap::pheatmap(fpkm_heatmap_tb,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   fontsize_row = 10, 
                   fontsize_col = 30,
                   border_color = NA,
                   angle_col = 45,
                   # col = colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))(20),
                   col = viridis(20),
                   # color = c("white", "white", "white", "white", str_c("grey", 99:1), "black", "black"),
                   # breaks = c(-1, 0, seq(1, 5.5, length = 100), 6),
                   # cellwidth = 50, 
                   # cellheight = 50, 
                   height = 50,
                   width = 50,
                   gaps_col = which.gap.col, 
                   filename = file.path(outpath, str_c("heatmap.developmental_profile.ensembl", ensembl_version, "KRAB_domain.pfam01352.pdf", sep = ".")))


### plot interactive heatmap using d3heatmap
# interactive heatmap
interactive_heatmap <- d3heatmap::d3heatmap(fpkm_heatmap_tb, 
                                            scale = "none", 
                                            Rowv = T, 
                                            Colv = F, 
                                            color = viridis(20), 
                                            cexCol = 0.2, 
                                            cexRow = 1, 
                                            show_grid = T, 
                                            anim_duration = 0, 
                                            cellnote = round(((10^fpkm_heatmap_tb) - 0.001), 2))

# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(interactive_heatmap),
                        file = file.path(outpath, str_c("heatmap.developmental_profile.ensembl", ensembl_version, "KRAB_domain.pfam01352.clustered.html", sep = ".")),
                        selfcontained = T)
