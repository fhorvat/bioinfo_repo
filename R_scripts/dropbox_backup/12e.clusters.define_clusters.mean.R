### INFO: 
### DATE: Fri Aug 31 15:10:23 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/smallRNA_clusters")

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
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### for each region sum weighted counts of all reads
GetOverlaps <- function(reg1, reg2, colname = NULL){
  
  if(is.null(colname)){
    stop("Colname needs to be defined")
  }
  
  fo <- data.table(as.matrix(findOverlaps(reg1, reg2)), ignore.strand = T)
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
l.files <- list.files(inpath, full.names = T, pattern = ".*\\.weighted_reads\\.RDS")

# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/library_size.Eliska_mESC_MosIR.csv"

# sample table path 
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Documentation/Eliska_mESC_MosIR.sampleTable.csv"

# annotation RDS path
annot_path <- file.path(inpath, "l.annot.RDS")

######################################################## READ DATA
# read bams from count matrices
bams <- 
  purrr::map(l.files, readRDS) %>% 
  magrittr::set_names(., l.files %>% basename(.) %>% str_remove(., "\\.SE\\.mis_0\\.18to30nt\\.weighted_reads\\.RDS"))

# read library size df
library_size_df <- readr::read_csv(file = library_size_path) 

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read annotation
l.annot <- readRDS(annot_path) 

######################################################## MAIN CODE
#### prepare tables ####
# filter sample table
sample_table %<>% 
  dplyr::select(sample_id, genotype, transfection) %>% 
  tidyr::unite(gen_trans, genotype, transfection, sep = ".") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE"))

# set library size to named vector
library_size <- 
  library_size_df$library_size.mis_0.18to30nt %>% 
  magrittr::set_names(., library_size_df$sample_id)

#### find clusters in bam files ####
# set cutoff values
rpkm1 <- 3
width <- 50
clust_list <- purrr::map(names(bams), function(x) FindClusters(reads = bams[[x]], total = library_size[[x]], rpkm1 = rpkm1, clust.width = width, e = 1e6))
names(clust_list) <- names(bams)


#### FIND AND PLOT CLUSTERS #### 
## X-axis = RS10_Mos or RSP_Mos, Y-axis = RS7_Mos

#### merge replicate clusters (intersect) ####
## loop through X-axis samples
for(x in c("RS10_Mos", "RSP_Mos")){
  
  # get names of samples
  x_samples <- names(bams)[str_detect(names(bams), x)]
  y_samples <- names(bams)[str_detect(names(bams), "RS7_Mos")]
  
  # merge all replicate
  x_clust <- purrr::reduce(clust_list[x_samples], ClusterMerge, method = "intersect")
  y_clust <- purrr::reduce(clust_list[y_samples], ClusterMerge, method = "intersect")
  
  
  #### merge all clusters (union) ####
  # merge all
  clust_all <- ClusterMerge(x_clust, y_clust, method = "union")
  
  # get clusters in table, annotate
  clust_all_df <- 
    clust_all %>% 
    as.data.frame(.) %>% 
    as.tibble(.) %>% 
    tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
    dplyr::mutate(class = AnnotateRegions(clust_all, l.annot))
  
  
  #### count reads in clusters ####
  # for each bam count reads in clusters
  d.cnts <- 
    do.call(cbind, lapply(bams, function(x) GetOverlaps(clust_all, x, "nh"))) %>% 
    as.tibble(.)
  
  # normalize counts for library size
  d.cnts.norm <- 
    d.cnts %>% 
    dplyr::mutate(cluster_id = 1:nrow(.)) %>% 
    tidyr::gather(key = "sample_id", value = count, -cluster_id) %>% 
    dplyr::left_join(., library_size_df, by = "sample_id") %>% 
    dplyr::mutate(library_size = library_size.mis_0.18to30nt / 1e6,
                  rpm = round((count / library_size), 3)) %>% 
    dplyr::select(cluster_id, sample_id, rpm) %>% 
    tidyr::spread(sample_id, rpm) %>% 
    dplyr::select(-cluster_id)
  
  
  #### output clusters table ####
  # for each loci get percentage of reads mapping to +/- strand
  cp <- GetOverlaps(clust_all, bams[[x_samples[1]]][strand(bams[[x_samples[1]]]) == "+"], "nh")
  cm <- GetOverlaps(clust_all, bams[[x_samples[1]]][strand(bams[[x_samples[1]]]) == "-"], "nh")
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
    readr::write_csv(., file.path(outpath, str_c("clusters.Eliska_mESC_MosIR", str_c(str_remove(x, "_Mos"), "_RS7"), "union", "rpm", rpkm1, "width", width, "csv", sep = ".")))
  
  # get long data.frame
  do_long <- 
    do %>% 
    dplyr::select(-c(width, strand, class)) %>% 
    tidyr::gather(sample_id, rpm, -coordinates) %>% 
    dplyr::mutate(sample = str_remove_all(sample_id, "^s_|_r[0-9]+$")) %>% 
    dplyr::group_by(coordinates, sample) %>% 
    dplyr::summarise(rpm_mean = round(mean(rpm), 3)) %>% 
    dplyr::ungroup(.) %>% 
    tidyr::gather(summary, rpm, rpm_mean) %>% 
    tidyr::unite(sample, sample, summary, sep = ".") %>% 
    tidyr::spread(sample, rpm) %>% 
    dplyr::left_join(., clust_all_df, by = "coordinates") %>% 
    dplyr::select(coordinates, width:class, everything())
  
  # get and save mean values across genotype
  do_mean <- 
    do_long %>% 
    dplyr::select_at(vars(coordinates, class, contains("mean"))) %T>% 
    readr::write_csv(., file.path(outpath, str_c("clusters.Eliska_mESC_MosIR", str_c(str_remove(x, "_Mos"), "_RS7"), "union", "rpm", rpkm1, "width", width, "mean", "csv", sep = ".")))
  
  # get cluster coordinates
  clusters_gr <- 
    do_long %>% 
    dplyr::select(coordinates, strand, class) %>% 
    tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
    GenomicRanges::GRanges(.)
  
  
  #### overlap with Dicer = ENSMUSG00000041415 ####
  # get Dicer1 coordinates
  dicer_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000041415")]
  
  # get hits overlaping with Dicer
  dicer_hits <- 
    findOverlaps(clusters_gr, dicer_gr) %>% 
    queryHits(.)
  
  
  #### overlap with Mos = ENSMUSG00000078365 ####
  # get Mos genomic coordinates
  mos_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000078365")]
  
  # get hits overlaping with genomic Mos
  mos_hits <-
    findOverlaps(clusters_gr, mos_gr) %>%
    queryHits(.)
  
  
  #### plot clusters ####
  # what samples to plot
  y <- "RS7_Mos"
  
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
    ggsave(filename = file.path(outpath, str_c("clusters.Eliska_mESC_MosIR", str_c(str_remove(x, "_Mos"), "_RS7"), "union", "rpm", rpkm1, 
                                               "width", width, "meanRPM", filter_name, "png", sep = ".")), 
           width = 12, height = 12)
  
}

