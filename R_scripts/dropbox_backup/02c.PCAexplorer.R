library("pcaExplorer")
install.packages("rlang")

se <- readRDS(file = file.path(outpath, "GRCm38.89.reducedExons.summarizedOverlaps.RDS")) 
sample_table_df <- 
  sample_table %>% 
  as.data.frame(.) %>% 
  magrittr::set_rownames(., .$sample_id)

colData(se) <- DataFrame(sample_table_df)
dds_test <- DESeqDataSet(se, design= ~stage)
rld_test <- rlogTransformation(dds_test)

### PCA 
pcaplot(rld_test, intgroup = "stage", ntop = 1000,
        pcX = 1, pcY = 2, title = "mESC_oocytes_2018 dataset PCA - PC1 vs PC2") + 
  ggsave(filename = file.path(outpath, "results", str_c("test.png")), width = 10, height = 10)

### expression distribution
distro_expr(rld_test, plot_type = "density") +
  ggsave(filename = file.path(outpath, "results", str_c("test2.png")), width = 10, height = 10)

distro_expr(rld_test, plot_type = "violin") + 
  ggsave(filename = file.path(outpath, "results", str_c("test3.png")), width = 10, height = 10)

distro_expr(rld_test, plot_type = "boxplot") + 
  ggsave(filename = file.path(outpath, "results", str_c("test4.png")), width = 10, height = 10)


### genespca
# annotation test
anno_df_biomart <- get_annotation(dds = dds_test,
                                  biomart_dataset = "mmusculus_gene_ensembl",
                                  idtype = "ensembl_gene_id")

groups_test <- colData(dds_test)$stage
cols_test <- scales::hue_pal()(3)[groups_test]

# with many genes, do not plot the labels of the genes
genespca(rld_test,
         ntop = 5000,
         choices = c(1, 2),
         arrowColors = cols_test,
         groupNames = groups_test,
         alpha = 0.2,
         useRownamesAsLabels = FALSE,
         varname.size = 5) + 
  ggsave(filename = file.path(outpath, "results", str_c("test5.png")), width = 10, height = 10)

# with a smaller number of genes, plot gene names included in the annotation
genespca(rld_test,
         ntop = 100,
         choices = c(1,2),
         arrowColors = cols_test,
         groupNames = groups_test,
         alpha = 0.7,
         varname.size = 5,
         annotation = anno_df_biomart) +
  ggsave(filename = file.path(outpath, "results", str_c("test6.png")), width = 10, height = 10)
