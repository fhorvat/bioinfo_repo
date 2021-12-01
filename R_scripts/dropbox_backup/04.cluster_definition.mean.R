b### INFO: 
### DATE: Tue Jan 07 15:36:14 2020
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/mESC_MosIR/datasets/Jul_2018/Analysis/expression/smallRNA_clusters")

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
l.files <- list.files(inpath, full.names = T, pattern = ".*filt.RData")

# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/mESC_MosIR/datasets/Jul_2018/Data/Mapped/STAR_mm10/4_library_size/library_sizes.txt"

# sample table path 
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/mESC_MosIR/datasets/Jul_2018/Data/Documentation/Eliska_mESC_MosIR.sampleTable.csv"

# annotation RDS path
annot_path <- file.path(inpath, "l.annot.RDS")

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
library_size_df <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read annotation
l.annot <- readRDS(annot_path)

######################################################## MAIN CODE
### prepare tables
# tidy sample table
sample_table %<>% 
  dplyr::select(sample_id, sample_name)

# tidy library sizes 
library_size_df %<>% 
  dplyr::filter(str_detect(sample_id, pattern = "\\.21to23nt$")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.21to23nt$"))

# set library size to named vector
library_size <- 
  library_size_df$library_size %>% 
  magrittr::set_names(., library_size_df$sample_id)


### find clusters in bam files
# set cutoff values
rpkm1 <- 3
width <- 50
clust_list <- purrr::map(names(bams), function(x) FindClusters(bams[[x]], library_size[[x]], rpkm1, width, 1e6))
names(clust_list) <- names(bams)


### merge clusters in replicates - RS7 vs. RS10 with Mos IR
# merge clusters from RS7 replicates (intersect)
RS7_samples <- names(bams)[str_detect(names(bams), "RS7_Mos")]
RS7_clust <- 
  ClusterMerge(clust_list[[RS7_samples[1]]], clust_list[[RS7_samples[2]]], method = "intersect") %>% 
  ClusterMerge(., clust_list[[RS7_samples[3]]], method = "intersect")

# merge clusters from RS10 replicates (intersect)
RS10_samples <- names(bams)[str_detect(names(bams), "RS10_Mos")]
RS10_clust <- 
  ClusterMerge(clust_list[[RS10_samples[1]]], clust_list[[RS10_samples[2]]], method = "intersect") %>% 
  ClusterMerge(., clust_list[[RS10_samples[3]]], method = "intersect")

# union of all samples
clust_all <- ClusterMerge(RS7_clust, RS10_clust, method = "union")

# get clusters data.frame
clust_all_df <- 
  clust_all %>% 
  as.data.frame(.) %>% 
  as_tibble(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::mutate(class = AnnotateRegions(clust_all, l.annot))


### count reads in clusters
# for each bam count reads in clusters
d.cnts <- 
  do.call(cbind, lapply(bams, function(x) GetOverlaps(clust_all, x, "nh"))) %>% 
  as_tibble(.)

# normalize counts for library size
d.cnts.norm <- 
  d.cnts %>% 
  dplyr::mutate(cluster_id = 1:nrow(.)) %>% 
  tidyr::pivot_longer(cols = -cluster_id, names_to = "sample_id", values_to = "count") %>% 
  dplyr::left_join(., library_size_df, by = "sample_id") %>% 
  dplyr::mutate(library_size = round((library_size / 1e6), 3),
                rpm = count / library_size) %>% 
  dplyr::select(cluster_id, sample_id, rpm) %>% 
  tidyr::pivot_wider(id_cols = cluster_id, names_from = sample_id, values_from = rpm) %>% 
  dplyr::select(-cluster_id)


### output RPM tables
# for each loci get percentage of reads mapping to +/- strand
cp <- GetOverlaps(clust_all, bams[[RS7_samples[1]]][strand(bams[[RS7_samples[1]]]) == "+"], "nh")
cm <- GetOverlaps(clust_all, bams[[RS7_samples[1]]][strand(bams[[RS7_samples[1]]]) == "-"], "nh")
perc <- cp / (cp + cm)

# add strand to clusters
clust_all_df %<>% 
  dplyr::mutate(strand = "*", 
                strand = replace(strand, perc > 0.8, "+"), 
                strand = replace(strand, perc < 0.2, "-"))

# create data.frame, set strand, save
do <- 
  cbind(clust_all_df, d.cnts.norm) %>% 
  tibble::as_tibble(.) %T>%
  readr::write_csv(., file.path(outpath, str_c("clusters.mESC_MosIR", "RS7", "RS10", "union", "rpm", rpkm1, "width", width, "csv", sep = ".")))

# calculate mean RPMs
do_long <- 
  do %>% 
  dplyr::select(-c(width, strand, class)) %>% 
  tidyr::pivot_longer(cols = -coordinates, names_to = "sample_id", values_to = "rpm") %>% 
  dplyr::left_join(., sample_table, by = "sample_id") %>% 
  dplyr::group_by(coordinates, sample_name) %>% 
  dplyr::summarise(rpm_mean = round(mean(rpm), 3)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = coordinates, names_from = sample_name, values_from = rpm_mean) %>% 
  dplyr::left_join(., clust_all_df, by = "coordinates") %>% 
  dplyr::select(coordinates, width:class, everything())

# get and save mean values across genotype
do_mean <- 
  do_long %>% 
  dplyr::select(-width) %T>% 
  readr::write_csv(., file.path(outpath, str_c("clusters.mESC_MosIR", "RS7", "RS10", "union", "rpm", rpkm1, 
                                               "width", width, "mean_RPM", "csv", sep = ".")))


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


### plot clusters as MA plot
# what samples to plot
x_sample <- "RS10_Mos"
y_sample <- "RS7_Mos"

# table for plot
plot_df <- 
  do_long %>% 
  dplyr::mutate(class = ifelse((1:nrow(.) %in% mos_hits), "Mos", class)) %>% 
  dplyr::mutate(class = factor(class, levels = c("Mos", "miRNA", "TE", "mRNA", "RNA_other", "other"))) %>% 
  dplyr::filter(!(1:nrow(.) %in% dicer_hits)) %>%
  dplyr::filter(!(str_detect(coordinates, "pCAG-EGFP_MosIR"))) %>%
  dplyr::arrange(desc(class)) %>% 
  dplyr::select(x = x_sample, y = y_sample, class) %>%
  dplyr::mutate(mean = ((x + y) / 2),
                lfc = log2(x + 1) - log2(y + 1))

# set filter name
filter_name <- "filtered.MosIR.Dicer"

### plot png 
# create ggplot object
clusters_plot <- 
  ggplot() +
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = class), size = 5, shape = 20) +
  scale_x_log10(limits = c(1, 1e6), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-8, 8)) +
  scale_color_manual(values = c("Mos" = "red", "miRNA" = "cornflowerblue", "TE" = "gray20",
                                "mRNA" = "red", "RNA_other" = "gray60",
                                "other" = "gray80")) +
  guides(color = FALSE) +
  # xlab("average expression (RPM)") +
  # ylab("log2(RS10 / RS7_Mos)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(filename = file.path(outpath, str_c("clusters", "MA_plot", "mESC_MosIR", "RS7", "RS10", "union", "rpm", rpkm1, 
                                           "width", width, "mean_RPM", filter_name, "png", sep = ".")), 
       plot = clusters_plot,
       width = 20, height = 10)


## add label to Mos
# filter Mos
label_df <- 
  plot_df %>% 
  dplyr::filter(class == "Mos")

# add label to plot
clusters_plot_labeled <- 
  clusters_plot +
  geom_point(data = label_df, aes(x = mean, y = lfc), color = "red", size = 8, shape = 17) + 
  guides(color = "legend") +
  xlab("average expression (RPM)") +
  ylab("log2(RS10 / RS7_Mos)") + 
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13, angle = 90))

# save
ggsave(filename = file.path(outpath, str_c("clusters", "MA_plot", "mESC_MosIR", "RS7", "RS10", "union", "rpm", rpkm1, 
                                           "width", width, "mean_RPM", filter_name, "labeled", "png", sep = ".")),
       plot = clusters_plot_labeled,
       width = 20, height = 10)


# ### plot clusters as scatter plot
# # what samples to plot
# x_sample <- "RS10_Mos"
# y_sample <- "RS7_Mos"
# 
# # table for plot
# plot_df <- 
#   do_long %>% 
#   dplyr::mutate(class = ifelse((1:nrow(.) %in% mos_hits), "Mos", class)) %>% 
#   dplyr::mutate(class = factor(class, levels = c("Mos", "miRNA", "TE", "mRNA", "RNA_other", "other"))) %>% 
#   dplyr::filter(!(1:nrow(.) %in% dicer_hits)) %>%
#   dplyr::filter(!(str_detect(coordinates, "pCAG-EGFP_MosIR"))) %>%
#   dplyr::arrange(desc(class)) %>% 
#   dplyr::select(x = x_sample, y = y_sample, class) %>%
#   dplyr::mutate(x = log10(x + 1), 
#                 y = log10(y + 1))
# 
# # set filter name
# filter_name <- "filtered.MosIR.Dicer"
# 
# ### plot png 
# # create ggplot object
# clusters_plot <- 
#   ggplot() +
#   geom_point(data = plot_df, aes(x = x, y = y, color = class), size = 2) +
#   scale_color_manual(values = c("Mos" = "red", "miRNA" = "cornflowerblue", "TE" = "gray20",
#                                 "mRNA" = "red", "RNA_other" = "gray60",
#                                 "other" = "gray80")) +
#   scale_x_continuous(limits = c(0, 6)) +
#   scale_y_continuous(limits = c(0, 6)) + 
#   xlab(x_sample) +
#   ylab(y_sample) +
#   coord_fixed() + 
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # save
# ggsave(filename = file.path(outpath, str_c("clusters", "scatter_plot", "mESC_MosIR", "RS7", "RS10", "union", "rpm", rpkm1, 
#                                            "width", width, "mean_RPM", filter_name, "png", sep = ".")), 
#        plot = clusters_plot,
#        width = 12, height = 12)


