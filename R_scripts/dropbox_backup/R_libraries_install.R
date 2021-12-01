#### CRAN ####
# packages
cran_packages <- c("dplyr", "readr", "magrittr", "stringr",  "tibble", "ggplot2", 
                   "data.table", "tidyr", "readxl", "Hmisc", "glue", "janitor", "scales", "devtools", "ggrepel", 
                   "plotly", "VennDiagram", "pheatmap", "openxlsx")

# install packages
install.packages(cran_packages, repos = "https://cloud.r-project.org/")

# load packages 
lapply(cran_packages, library, character.only = TRUE)


#### BIOCONDUCTOR ####
# packages
bioconductor_packages <- c("GenomicAlignments", "GenomicFeatures", "GenomicRanges", "ShortRead", "pathview", "org.Hs.eg.db",
                           "org.Mm.eg.db", "Biostrings", "BiocParallel", "SummarizedExperiment", "Rsamtools", "rtracklayer", "GenomeInfoDb", 
                           "AnnotationDbi", "Biobase", 
                           # "BSgenome.Mmusculus.UCSC.mm10", 
                           "DESeq2", "GeneOverlap", 
                           "oligo", "affy", "gcrma", "mouse4302.db", "pd.mogene.2.0.st")

# install Bioconductor package (only first time)
install.packages("BiocManager")

# install packages
BiocManager::install(bioconductor_packages)

# load packages
lapply(bioconductor_packages, require, character.only = TRUE)


