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

library(Cairo)
library(data.table)
library(GenomicRanges)
library(plyr)

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
  ind <- c(names(l.annot), "Other")[mi]
  
  # return
  return(ind)
  
}


######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get gene info path
info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# get miRBase gff path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# list count matrices RData files
counts_path <- file.path(inpath, "filtered_bams")
l.files <- list.files(counts_path, full.names = T, pattern = ".*filt.RData")

# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/MosIR_expression/mosIR.T3T_DcrTrans_2011.counts_summary.xlsx"

# sample table path 
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Documentation/T3T_DcrTrans_2011.sample_table.csv"

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_path) 

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
library_size_df <- 
  xlsx::read.xlsx(file = library_size_path, sheetName = "library_sizes") %>% 
  tibble::as.tibble(.) %>% 
  tidyr::gather(library_type, library_size, -sample_id) %>% 
  tidyr::separate(library_type, into = c("mis", "library_type"), sep = "\\.") %>% 
  tidyr::unite(sample_id, sample_id, mis, sep = ".") %>% 
  tidyr::spread(library_type, library_size)

# read sample table
sample_table <- 
  readr::read_csv(sample_table_path) %>% 
  dplyr::select(sample_id, genotype, transfection) %>% 
  tidyr::unite(gen_trans, genotype, transfection, sep = ".")

######################################################## MAIN CODE
### prepare annotations - repeatMasker, exons of genes from ENSEMBL annotation, miRBase
# prepare exons
exons_gr %<>% 
  as.data.frame(.) %>% 
  tibble::as.tibble(.) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>%
  dplyr::filter(!str_detect(gene_biotype, "pseudogene|miRNA")) %>% 
  GenomicRanges::GRanges(.)

# prepare repeats
rmsk_gr <- 
  rmsk_df %>% 
  GenomicRanges::GRanges(.) %>% 
  reduce(.)

# prepare mature miRNA
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name")]
names(mcols(mirna_gr)) <- "gene_id"
mcols(mirna_gr)$gene_biotype <- "miRNA.mature"

# add annotation to list
l.annot <- list(miRNA = mirna_gr,
                Genes = exons_gr,
                Repeats = rmsk_gr)

# clean memory
# rm(mirna_gr, exons_gr, rmsk_gr, rmsk_df); gc()

### find clusters in bam files
# DicerO samples
dicero_bams <- c("s_T3T_DcrO_tran_r1.SE.mis_0", "s_T3T_WT_utran_r1.SE.mis_0")
names(dicero_bams) <- c("DicerO_r1", "DicerWT_r1")

# set library size to named vector
library_size <- 
  library_size_df %>% 
  dplyr::filter(sample_id %in% dicero_bams) %>% 
  dplyr::select(sample_id, library_size = `21to23nt_reads`)

library_size <- 
  library_size$library_size %>% 
  magrittr::set_names(., library_size$sample_id)

# set cutoff values
rpkm1 <- 3
width <- 50
dicero_clust_list <- purrr::map(dicero_bams, function(x) FindClusters(bams[[x]], library_size[[x]], rpkm1, width, 1e6))

# merge clusters from replicates (method either intersect or union)
dicero_clust <- dicero_clust_list[[1]]

# get clusters data.frame
dicero_clust_df <- 
  dicero_clust %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::mutate(class = AnnotateRegions(dicero_clust, l.annot))


### count reads in clusters
# for each bam count reads in clusters
d.cnts <- 
  do.call(cbind, lapply(bams, function(x) GetOverlaps(dicero_clust, x, "nh"))) %>% 
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
cp <- GetOverlaps(dicero_clust, bams[[dicero_bams[1]]][strand(bams[[dicero_bams[1]]]) == "+"], "nh")
cm <- GetOverlaps(dicero_clust, bams[[dicero_bams[1]]][strand(bams[[dicero_bams[1]]]) == "-"], "nh")
perc <- cp / (cp + cm)

# add strand to clusters
dicero_clust_df %<>% 
  dplyr::mutate(strand = "*", 
                strand = replace(strand, perc > 0.8, "+"), 
                strand = replace(strand, perc < 0.2, "-"))

# create data.frame, set strand, save
do <- 
  cbind(dicero_clust_df, d.cnts.norm) %>% 
  tibble::as.tibble(.)


## overlap with Dicer = ENSMUSG00000041415
# get Dicer1 coordinates
dicer_gr <- exons_gr[exons_gr$gene_id == "ENSMUSG00000041415"]

# get cluster coordinates
clusters_gr <- 
  do_long %>% 
  dplyr::select(coordinates, strand, class) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# get hits overlaping with Dicer
dicer_hits <- 
  findOverlaps(clusters_gr, dicer_gr) %>% 
  queryHits(.)


### plot clusters
# what samples to plot
# what samples to plot
x <- "s_ES_Dcr_Rsc_tran_r2.SE.mis_0"
y_list <- c("s_ES_WT_tran_r1.SE.mis_0", "s_ES_DcrSom_Rsc_tran_r1.SE.mis_0", "s_ES_DcrSom_Rsc_tran_r2.SE.mis_0")

for(y in y_list){
  
  # get names
  name_x <- names(dicero_bams[dicero_bams == x])
  name_y <- names(dicero_bams[dicero_bams == y])
  
  # table for plot
  plot_df <- 
    do %>% 
    dplyr::filter(!(1:nrow(.) %in% dicer_hits)) %>%
    dplyr::filter_at(vars(matches(str_c(x, "|", y))), any_vars(. > 50)) %>%
    dplyr::select(x = x, y = y, class) %>%
    dplyr::mutate(x = log10(x + 1), 
                  y = log10(y + 1)) 
  
  # set filter name
  filter_name <- "filtered.any50.Dicer"
  
  # plot and save
  ggplot() +
    geom_point(data = plot_df, aes(x = x, y = y, color = class), size = 3) +
    scale_color_manual(values = c("Genes" = "red", "miRNA" = "cornflowerblue", "Other" = "black", "Repeats" = "darkgray", "MosIR" = "green")) +
    scale_x_continuous(limits = c(0, 6)) +
    scale_y_continuous(limits = c(0, 6)) + 
    xlab(name_x) +
    ylab(name_y) +
    guides(color = FALSE) + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(outpath, str_c("clusters", "DicerO3", "rpm", rpkm1, "width", width, 
                                               name_x, "vs", name_y, "RPM", filter_name, "png", sep = ".")), 
           width = 10, height = 10)
  
}
