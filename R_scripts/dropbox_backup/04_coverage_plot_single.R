positionCoverageDFSingle <- function(coverage_counts, stage){
  # creating position matrix with all counts
  all_positions <- unique(unlist(lapply(coverage_counts, names)))
  position_matrix <- matrix(NA, 
                            nrow = length(all_positions), 
                            ncol = length(coverage_counts), 
                            dimnames = list(all_positions, c(1:length(coverage_counts))))
  for (i in seq_along(coverage_counts)) {
    position_matrix[names(coverage_counts[[i]]), i] <- coverage_counts[[i]]
  }
  position_matrix <- position_matrix[complete.cases(position_matrix), ]
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000"):which(rownames(position_matrix) == "150000"), ]
  return(position_matrix)
}

positionCoverageDFOneElementAllStages <- function(x){
  one_element <- all_coverage_counts_df[, x, drop = F]
  one_element$pos <- rep(seq(-150000, 150000), 5)
  one_element$stage <- rep(all_stages_names, each = 300001)
  one_element$stage <- factor(one_element$stage, levels = c("GV", "1C", "2C", "2C_aphi", "4C"))
  one_element$element <- ifelse((one_element$pos >= 0 & one_element$pos <= width(all_elements_original_gr[all_elements_original_gr$fullName == colnames(one_element)[1]])), 
                                "in_element", 
                                "out_element")
  colnames(one_element)[1] <- "count"
  one_element$count[one_element$count == 0] <- NA
  one_element$pos[is.na(one_element$count)] <- NA
  return(one_element)
}

singleElementPlotByStages <- function(x){
  plot_df <- all_elements_coverage_list[[x]]
  single_element_facet_ggplot <- 
    ggplot(plot_df, aes(x = pos, y = count)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
    scale_x_continuous(limits = c(-150000, 150000)) + 
    scale_fill_manual(values = c("red", "black")) +
    facet_grid(stage ~ .) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_blank()) +
    coord_cartesian(ylim = c(0, 50)) +
    ggtitle(names(all_elements_coverage_list[x])) 
}

#########################################################################################################
all_stages_counts <- list(all_coverage_GV, 
                          all_coverage_1C, 
                          all_coverage_2C, 
                          all_coverage_2C_aphi, 
                          all_coverage_4C)
all_stages_names <- c("GV", "1C", "2C", "2C_aphi", "4C")
all_coverage_counts <- lapply(X = 1:length(all_stages_counts), 
                              function(X) positionCoverageDFSingle(coverage_counts = all_stages_counts[[X]], stage = all_stages_names[X]))

all_coverage_counts_df <- as.data.frame(do.call(rbind, all_coverage_counts))
rownames(all_coverage_counts_df) <- NULL
colnames(all_coverage_counts_df) <- names(all_coverage_GV)

#########################################################################################################
rm(list = ls()[grepl("all_coverage_counts$", ls())])

# list of elements with data.frames by stages (long format)
all_elements_coverage_list <- lapply(1:ncol(all_coverage_counts_df), positionCoverageDFOneElementAllStages)

# original data ranges ordered by names of all_elements_coverage_list
all_elements_original_ordered <- left_join(data.frame(fullName = names(all_coverage_GV)), all_elements_original,by = "fullName")

names(all_elements_coverage_list) <- paste0(all_elements_original_ordered$seqnames, "\t",
                                            all_elements_original_ordered$start + 1.5e5, "\t", 
                                            all_elements_original_ordered$end - 1.5e5, "\t strand: ", 
                                            all_elements_original_ordered$strand, "\t element: ", 
                                            all_elements_original_ordered$repName, "\n",
                                            names(all_coverage_GV))

#########################################################################################################
# plot MT2 solo
MT2_solo_plot_list <- lapply(grep("MT2_Mm", names(all_elements_coverage_list)), singleElementPlotByStages) 
pdf("MT2_Mm_top100_ylim_50_fpkm_1_filtered_single_coverage.pdf", width = 30, height = 10)
invisible(lapply(MT2_solo_plot_list, print))
dev.off()

# plot ORR1A0 solo
ORR1A0_plot_list <- lapply(grep("ORR1A0", names(all_elements_coverage_list)), singleElementPlotByStages)
pdf("ORR1A0_top100_ylim_50_fpkm_1_filtered_single_coverage.pdf", width = 30, height = 10)
invisible(lapply(ORR1A0_plot_list, print))
dev.off()

# plot MT2 full
MT2_full_plot_list <- lapply(grep("MT2_full", names(all_elements_coverage_list)), singleElementPlotByStages) 
pdf("MT2_full_top100_ylim_50_fpkm_1_filtered_single_coverage.pdf", width = 30, height = 10)
invisible(lapply(MT2_full_plot_list, print))
dev.off()

rm(list = ls()[grepl("all_elements_coverage_list|_plot_list", ls())])


