### INFO: 
### DATE: Fri Aug 31 15:10:23 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/T3T_DcrTrans_2011")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(data.table)
library(GenomicRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
### assign RData file to variable
Assigner <- function(`.path`, `.name`){
  
  if(!is.character(`.path`) | !is.character(`.name`)){
    stop('Both arguments should be characters!')
  }
  
  load(`.path`)
  
  assign(`.name`, get(ls()[1]), parent.frame())
  
}


### for each region sum weighted counts of all reads
GetOverlaps <- function(reg1, reg2, colname = NULL){
  
  if(is.null(colname)){
    stop("Colname needs to be defined")
  }
  
  fo <- data.table(as.matrix(findOverlaps(reg1, reg2)), ignore.strand=T)
  fo$weight <- values(reg2)[[colname]][fo$subjectHits]
  fo <- fo[, sum(weight), by = queryHits]
  v <- rep(0, length(reg1))
  v[fo$queryHits] <- fo$V1
  
  return(v)
  
}


### finds the clusters in short RNA data
FindClusters <- function(reads, total, rpkm1 = 5, clust.width = 50, e = 1e6){
  
  # reduce regions
  dregs <- reduce(reads, ignore.strand = T)
  
  # get weighted counts in each region
  co <- GetOverlaps(dregs, reads, "nh")
  
  # filtering based on RPKM
  rpkm <- (co) * (e / total)
  rind <- rpkm > rpkm1
  wregs <- dregs[rind]
  values(wregs)$counts <- co[rind]
  values(wregs)$rpkm <- rpkm[rind]
  
  # join cluster closer than clust.width
  a.regs <- reduce(resize(wregs, width = width(wregs) + clust.width, fix = "center"), ignore.strand = T)
  a.regs <- resize(a.regs, width = width(a.regs) - clust.width, fix = "center")
  values(a.regs)$counts <- GetOverlaps(a.regs, reads, "nh")
  values(a.regs)$rpkm.c <- (values(a.regs)$counts) * (e / total)
  regs.sel <- a.regs
  
  # return 
  return(reg.sel = regs.sel)
  
}


### merges sets of cluster
ClusterMerge <- function(regs1, regs2, method = "intersect", ignore.strand = T){
  
  cat("method:", method, "\n")
  
  # intersect
  if(method == "intersect"){
    
    fo <- as.matrix(findOverlaps(regs1, regs2))
    o <- reduce(c(regs1[fo[, 1]], regs2[fo[, 2]]), ignore.strand = ignore.strand)
    
  }
  
  # union
  if(method == "union"){
    
    o <- reduce(c(regs1, regs2), ignore.strand = ignore.strand)
    
  }
  
  # return
  return(regs.sel = o)
  
}


### hierarchically annotate regions with list of annotations  
AnnotateRegions <- function(regs, l.annot){
  
  # count overlaps between regions and annotation
  d <- do.call(cbind, lapply(l.annot, function(x) suppressWarnings(countOverlaps(regs, x))))
  
  # set all regions which overlap more than 1 annotation of same type to 1
  d[d > 1] <- 1
  
  # multiply numbers in each column with index of that column
  d <- t(t(d) * 1:ncol(d))
  
  # set all regions which don't overlap with any annotation to highest number
  d[d == 0] <- max(d) + 1
  
  
  # for each region get lowest numbered annotation
  da <- t(data.frame(d))
  mi <- do.call(pmin, split(da, 1:ncol(d)))
  
  # get names of annotation for each region
  ind <- c(names(l.annot), "not_annotated")[mi]
  
  # return
  return(ind)
  
}


######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list count matrices RData files
counts_path <- file.path(inpath, "filtered_bams")
l.files <- list.files(counts_path, full.names = T, pattern = ".*filt.RData")

# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/MosIR_expression/mosIR.T3T_DcrTrans_2011.counts_summary.xlsx"

# sample table path 
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Documentation/T3T_DcrTrans_2011.sample_table.csv"

# annotation RDS path
annot_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/l.annot.RDS"

######################################################## READ DATA
# read count matrices
l.mat <- list()
for(i in 1:length(l.files)){
  
  l.file <- l.files[i]
  name <- stringr::str_remove(basename(l.file), ".filt.RData")
  print(name)
  Assigner(l.file, "l")
  l.mat[[name]] <- l
  
}

# read bams from count matrices
bams <- lapply(l.mat, function(x) x$reads)

# read library size df
library_size_df <- xlsx::read.xlsx(file = library_size_path, sheetName = "library_sizes") 

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read annotation
l.annot <- readRDS(annot_path)

######################################################## MAIN CODE
### prepare tables
# filter sample table
sample_table %<>% 
  dplyr::select(sample_id, genotype, transfection) %>% 
  tidyr::unite(gen_trans, genotype, transfection, sep = ".")

# get library size from xlsx
library_size_df %<>% 
  tibble::as.tibble(.) %>% 
  tidyr::gather(library_type, library_size, -sample_id) %>% 
  tidyr::separate(library_type, into = c("mis", "library_type"), sep = "\\.") %>% 
  tidyr::unite(sample_id, sample_id, mis, sep = ".") %>% 
  tidyr::spread(library_type, library_size)

# set library size to named vector
library_size <- 
  library_size_df$`21to23nt_reads` %>% 
  magrittr::set_names(., library_size_df$sample_id)

### find clusters in bam files
# set cutoff values
rpkm1 <- 3
width <- 50
clust_list <- purrr::map(names(bams), function(x) FindClusters(bams[[x]], library_size[[x]], rpkm1, width, 1e6))
names(clust_list) <- names(bams)

# merge clusters from DicerO replicates (intersect)
dicerO_samples <- names(bams)[str_detect(names(bams), "DcrO")]
dicerO_clust <- 
  ClusterMerge(clust_list[[dicerO_samples[1]]], clust_list[[dicerO_samples[2]]], method = "intersect") %>% 
  ClusterMerge(., clust_list[[dicerO_samples[3]]], method = "intersect")

# merge clusters from DicerS replicates (intersect)
dicerS_samples <- names(bams)[str_detect(names(bams), "DcrS")]
dicerS_clust <- 
  ClusterMerge(clust_list[[dicerS_samples[1]]], clust_list[[dicerS_samples[2]]], method = "intersect") %>% 
  ClusterMerge(., clust_list[[dicerS_samples[3]]], method = "intersect")

# merge clusters from DicerS replicates (intersect)
dicerWT_samples <- names(bams)[str_detect(names(bams), "WT")]
dicerWT_clust <- 
  ClusterMerge(clust_list[[dicerWT_samples[1]]], clust_list[[dicerWT_samples[2]]], method = "intersect") %>% 
  ClusterMerge(., clust_list[[dicerWT_samples[3]]], method = "intersect")

# y-axis clusters
clusters_yaxis <- list(dicerS_clust, dicerWT_clust)
names(clusters_yaxis) <- c("DicerS", "DicerWT")


### get scatterplots of cluster RPM expression
for(clust_name in names(clusters_yaxis)){
  
  # get cluster for y-axis
  clust_y <- clusters_yaxis[[clust_name]]
  
  # union of all samples
  clust_all <- ClusterMerge(dicerO_clust, clust_y, method = "union")
  
  # get clusters data.frame
  clust_all_df <- 
    clust_all %>% 
    as.data.frame(.) %>% 
    as.tibble(.) %>% 
    tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
    dplyr::mutate(class = AnnotateRegions(clust_all, l.annot))
  
  ### count reads in clusters
  # for each bam count reads in clusters
  d.cnts <- 
    do.call(cbind, lapply(bams, function(x) GetOverlaps(clust_all, x, "nh"))) %>% 
    as.tibble(.)
  
  # normalize counts for library size
  d.cnts.norm <- 
    d.cnts %>% 
    dplyr::mutate(cluster_id = 1:nrow(.)) %>% 
    tidyr::gather(key = "sample_id", value = count, -cluster_id) %>% 
    dplyr::left_join(., library_size_df %>% dplyr::select(sample_id, library_size = `21to23nt_reads`), by = "sample_id") %>% 
    dplyr::mutate(library_size = round((library_size / 1e6), 4),
                  rpm = count / library_size) %>% 
    dplyr::select(cluster_id, sample_id, rpm) %>% 
    tidyr::spread(sample_id, rpm) %>% 
    dplyr::select(-cluster_id)
  
  ### output DicerO table
  # for each loci get percentage of reads mapping to +/- strand
  cp <- GetOverlaps(clust_all, bams[[dicerO_samples[1]]][strand(bams[[dicerO_samples[1]]]) == "+"], "nh")
  cm <- GetOverlaps(clust_all, bams[[dicerO_samples[1]]][strand(bams[[dicerO_samples[1]]]) == "-"], "nh")
  perc <- cp / (cp + cm)
  
  # add strand to clusters
  clust_all_df %<>% 
    dplyr::mutate(strand = "*", 
                  strand = replace(strand, perc > 0.8, "+"), 
                  strand = replace(strand, perc < 0.2, "-"))
  
  # create data.frame, set strand, save
  do <- 
    cbind(clust_all_df, d.cnts.norm) %>% 
    tibble::as.tibble(.) %T>%
    readr::write_csv(., file.path(outpath, str_c("clusters.3T3", "DicerO", clust_name, "union", "rpm", rpkm1, "width", width, "csv", sep = ".")))
  
  # get long data.frame
  do_long <- 
    do %>% 
    dplyr::select(-c(width, strand, class)) %>% 
    tidyr::gather(sample_id, rpm, -coordinates) %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, ".mis_.*")) %>% 
    dplyr::left_join(., sample_table, by = "sample_id") %>% 
    dplyr::group_by(coordinates, gen_trans) %>% 
    dplyr::summarise(rpm_mean = round(mean(rpm), 3)) %>% 
    dplyr::ungroup(.) %>% 
    tidyr::gather(summary, rpm, rpm_mean) %>% 
    tidyr::unite(gen_trans, gen_trans, summary, sep = ".") %>% 
    tidyr::spread(gen_trans, rpm) %>% 
    dplyr::left_join(., clust_all_df, by = "coordinates") %>% 
    dplyr::select(coordinates, width:class, everything())
  
  # get and save mean values across genotype
  do_mean <- 
    do_long %>% 
    dplyr::select_at(vars(coordinates, class, contains("mean"))) %T>% 
    readr::write_csv(., file.path(outpath, str_c("clusters.3T3", "DicerO", clust_name, "union", "rpm", rpkm1, "width", width, "mean", "csv", sep = ".")))
  
  ## overlap with Dicer = ENSMUSG00000041415
  # get Dicer1 coordinates
  dicer_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000041415")]
  
  # get cluster coordinates
  clusters_gr <- 
    do_long %>% 
    dplyr::select(coordinates, strand, class) %>% 
    tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
    GenomicRanges::GRanges(.)
  
  # get hits overlaping with Dicer
  dicer_hits <- 
    findOverlaps(clusters_gr, dicer_gr) %>% 
    queryHits(.)
  
  ## overlap with Mos = ENSMUSG00000078365
  # get Mos genomic coordinates
  mos_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000078365")]
  
  # get hits overlaping with genomic Mos
  mos_hits <-
    findOverlaps(clusters_gr, mos_gr) %>%
    queryHits(.)
  
  ### plot clusters
  # what samples to plot
  x <- "DcrO.tran"
  if(clust_name == "DicerS"){
    y <- "DcrS.tran"
  }else{
    y <- "WT.utran"
  }
  
  # prepare vectors to plot (DcrO.tran, DcrS.tran, WT.utran)
  mean.x <- str_c(x, ".rpm_mean")
  mean.y <- str_c(y, ".rpm_mean")
  
  # table for plot
  plot_df <- 
    do_long %>% 
    dplyr::mutate(class = ifelse((1:nrow(.) %in% mos_hits), "Mos", class)) %>% 
    dplyr::mutate(class = factor(class, levels = c("Mos", "miRNA", "TE", "mRNA", "RNA_other", "other"))) %>% 
    dplyr::filter(!(1:nrow(.) %in% dicer_hits)) %>%
    dplyr::filter(!(str_detect(coordinates, "pCAG-EGFP_MosIR"))) %>%
    dplyr::arrange(desc(class)) %>% 
    dplyr::select(x = mean.x, y = mean.y, class) %>%
    dplyr::mutate(x = log10(x + 1), 
                  y = log10(y + 1))
  
  # set filter name
  filter_name <- "filtered.MosIR.Dicer"
  
  ### plot png 
  ggplot() +
    geom_point(data = plot_df, aes(x = x, y = y, color = class), size = 2) +
    scale_color_manual(values = c("Mos" = "red", "miRNA" = "cornflowerblue", "TE" = "gray20",
                                  "mRNA" = "red", "RNA_other" = "gray60",
                                  "other" = "gray80")) +
    scale_x_continuous(limits = c(0, 6)) +
    scale_y_continuous(limits = c(0, 6)) + 
    xlab(x) +
    ylab(y) +
    coord_fixed() + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(outpath, str_c("clusters.3T3", "DicerO", clust_name, "union", "rpm", rpkm1, 
                                               "width", width, "meanRPM", filter_name, "png", sep = ".")), 
           width = 12, height = 12)
  
}


