### interactive 3D PCA plot
# calculates pca
pca <-
  rlog_df %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# data for plot
pca_3Dplot <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         PC3 = pca$x[, 3],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "), 
                stage = factor(stage, levels = c("GV", "MII", "1C")), 
                genotype = factor(genotype, levels = c("WT", "KO")))

# make interactive MA plot
p <-
  plotly::plot_ly(data = pca_3Dplot,
                  x = ~PC1,
                  y = ~PC2,
                  z = ~PC3,
                  symbol = ~stage,
                  symbols = c("square", "diamond", "circle"),
                  color = ~genotype,
                  colors = c("#1a75ff", "red3")) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = str_c("PC1: ", round(percentVar[1] * 100), "% variance"), range = c(-40, 40)),
                      yaxis = list(title = str_c("PC2: ", round(percentVar[2] * 100), "% variance"), range = c(-40, 40)), 
                      zaxis = list(title = str_c("PC3: ", round(percentVar[3] * 100), "% variance"), range = c(-40, 40))))


# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "results", str_c("PCAplot.3D.CNOT6L.GRCm38.89.PC1_2_3.rlog.html")),
                        selfcontained = T)