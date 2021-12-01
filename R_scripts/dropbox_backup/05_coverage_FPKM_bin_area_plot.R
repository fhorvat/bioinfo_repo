plotCumFPKMBinArea <- function(tile_width, all_coverage_counts_df_merged, which_element){
  
  element_width <- max(width(MT2_ORR1A0_original_gr[MT2_ORR1A0_original_gr$repName == which_element]))
  
  # upstream
  all_coverage_counts_df_merged_upstream <- all_coverage_counts_df_merged[all_coverage_counts_df_merged$pos >= -150000 & all_coverage_counts_df_merged$pos < 0, ]
  all_coverage_counts_df_merged_upstream <- 
    all_coverage_counts_df_merged_upstream %>%
    group_by(indx = gl(ceiling(nrow(all_coverage_counts_df_merged_upstream)/tile_width), tile_width, nrow(all_coverage_counts_df_merged_upstream))) %>%
    select(2:6) %>%
    summarise_each(funs(sum))
  all_coverage_counts_df_merged_upstream$position <- paste(seq(150000, tile_width, -tile_width), "upstream", sep = "_")
  
  # element
  all_coverage_counts_df_merged_element <- all_coverage_counts_df_merged[all_coverage_counts_df_merged$pos >= 0 & all_coverage_counts_df_merged$pos <= element_width, ]
  all_coverage_counts_df_merged_element <- 
    all_coverage_counts_df_merged_element %>%
    select(2:6) %>%
    summarise_each(funs(sum))
  all_coverage_counts_df_merged_element$position <- which_element 
  all_coverage_counts_df_merged_element$indx <- 1
  all_coverage_counts_df_merged_element <- all_coverage_counts_df_merged_element[, c(7, 1:6)]
  
  # downstream
  all_coverage_counts_df_merged_downstream <- all_coverage_counts_df_merged[all_coverage_counts_df_merged$pos > element_width & 
                                                                              all_coverage_counts_df_merged$pos <= (150000 + element_width), ]
  all_coverage_counts_df_merged_downstream <- 
    all_coverage_counts_df_merged_downstream %>%
    group_by(indx = gl(ceiling(nrow(all_coverage_counts_df_merged_downstream)/tile_width), tile_width, nrow(all_coverage_counts_df_merged_downstream))) %>%
    select(2:6) %>%
    summarise_each(funs(sum))
  all_coverage_counts_df_merged_downstream$position <- paste(seq(0, (150000 - tile_width), tile_width), "downstream", sep = "_")
  
  # creating one data.frame with mt2 and upstream/downstream positions
  all_coverage_counts_df_sum <- rbind(all_coverage_counts_df_merged_upstream, all_coverage_counts_df_merged_element, all_coverage_counts_df_merged_downstream)
  all_coverage_counts_df_sum <- all_coverage_counts_df_sum[, -1]
  
  # calculating FPKM for bins
  tile_kb <- tile_width / 1000
  all_coverage_counts_df_sum_fpkm <- all_coverage_counts_df_sum
  all_coverage_counts_df_sum_fpkm[, "GV"] <- all_coverage_counts_df_sum_fpkm[, "GV"] / (number_of_reads["GV"] * tile_kb)
  all_coverage_counts_df_sum_fpkm[, "1C"] <- all_coverage_counts_df_sum_fpkm[, "1C"] / (number_of_reads["1C"] * tile_kb)
  all_coverage_counts_df_sum_fpkm[, "2C"] <- all_coverage_counts_df_sum_fpkm[, "2C"]  / (number_of_reads["2C"] * tile_kb)
  all_coverage_counts_df_sum_fpkm[, "2C_aphi"] <- all_coverage_counts_df_sum_fpkm[, "2C_aphi"] / (number_of_reads["2C_aphi"] *  tile_kb)
  all_coverage_counts_df_sum_fpkm[, "4C"] <- all_coverage_counts_df_sum_fpkm[, "4C"] / (number_of_reads["4C"] * tile_kb)
  all_coverage_counts_df_sum_fpkm <- as.data.frame(all_coverage_counts_df_sum_fpkm)
  
  # melting data.frame to long format for ggplot
  all_coverage_counts_df_sum_fpkm_melt <- melt(all_coverage_counts_df_sum_fpkm, id.vars = "position")
  colnames(all_coverage_counts_df_sum_fpkm_melt) <- c("position", "stage", "fpkm")
  all_coverage_counts_df_sum_fpkm_melt$element <- ifelse((all_coverage_counts_df_sum_fpkm_melt$position == which_element),"in_element",  "out_element")
  all_coverage_counts_df_sum_fpkm_melt$width <- ifelse((all_coverage_counts_df_sum_fpkm_melt$position == which_element), element_width, tile_width)
  all_coverage_counts_df_sum_fpkm_melt$pos <- rep(c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, (150000 + element_width - tile_width), tile_width)), 5)
  all_coverage_counts_df_sum_fpkm_melt <- all_coverage_counts_df_sum_fpkm_melt[, c("pos", "fpkm", "mt2", "stage", "position", "width")]
  
  ## output .csv
  # all_coverage_counts_df_sum_fpkm_output <- all_coverage_counts_df_sum_fpkm
  # all_coverage_counts_df_sum_fpkm_output$pos <- c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, (150000 + element_width - tile_width), tile_width))
  # write.csv(all_coverage_counts_df_sum_fpkm_output, paste0("MT2_top100_AllStages_cummulativeFPKM_bins", tile_width, ".csv"), row.names = F)
  
  # order by stage
  stage_order <- rep(c("2C_aphi", "2C", "1C", "4C", "GV"), each = (nrow(all_coverage_counts_df_sum_fpkm_melt) / 5))
  all_coverage_counts_df_sum_fpkm_melt$stage <- factor(all_coverage_counts_df_sum_fpkm_melt$stage, levels = stage_order) 
  all_coverage_counts_df_sum_fpkm_melt <- all_coverage_counts_df_sum_fpkm_melt[order(all_coverage_counts_df_sum_fpkm_melt$stage), ]  
  
  # # bin plot
  all_coverage_counts_df_sum_fpkm_melt$stage <- factor(all_coverage_counts_df_sum_fpkm_melt$stage, levels = c("GV", "4C", "1C", "2C", "2C_aphi"))
  ggplot(all_coverage_counts_df_sum_fpkm_melt, aes(x = pos, y = fpkm)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + width, ymax = fpkm, fill = stage)) +
    scale_fill_manual(values = c("grey80", "grey60", "grey40", "grey20", "orange"), 
                      breaks = c("GV", "1C", "2C", "2C_aphi", "4C")) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width), 
                       breaks = c(seq(-150000, -10000, 10000), 0, seq(element_width, (140000 + element_width), 10000)), 
                       labels = c(seq(150, 10, -10), "MuERV", seq(0, 140, 10)), 
                       name = "bin") + 
    scale_y_continuous(name = "Cummulative FPKM") +  
    coord_cartesian(ylim = c(0, 2000)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = paste0(which_element, "_top100_AllStages_cummulativeFPKM_bins", tile_width, "_binPlot.pdf"), width = 30, height = 10)
  
  # area plot
  all_coverage_counts_df_sum_fpkm_melt$stage <- factor(all_coverage_counts_df_sum_fpkm_melt$stage, levels = rev(c("GV", "4C", "1C", "2C", "2C_aphi")))
  ggplot(all_coverage_counts_df_sum_fpkm_melt, aes(x = pos, y = fpkm)) +
    geom_area(aes(fill = stage), position = "identity") +
    scale_fill_manual(values = rev(c("grey80", "grey60", "grey40", "grey20", "orange")), 
                      breaks = c("GV", "1C", "2C", "2C_aphi", "4C")) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width), 
                       breaks = c(seq(-150000, -10000, 10000), 0, seq(element_width, (140000 + element_width), 10000)), 
                       labels = c(seq(150, 10, -10), "MuERV", seq(0, 140, 10)), 
                       name = "bin") + 
    scale_y_continuous(name = "Cummulative FPKM") + 
    coord_cartesian(ylim = c(0, 2000)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = paste0(which_element, "_top100_AllStages_cummulativeFPKM_bins", tile_width, "_areaPlot.pdf"), width = 30, height = 10)
}

# MT2 plot
MT2_all_coverage_counts_list_merged <- lapply(X = 1:length(all_stages_counts), FUN = function(X) MT2_all_coverage_counts[[X]][, 1:2])
MT2_all_coverage_counts_df_merged <- Reduce(function(...) merge(..., by = 'pos', all.x = TRUE), MT2_all_coverage_counts_list_merged)
colnames(MT2_all_coverage_counts_df_merged) <- c("pos", "GV", "1C", "2C", "2C_aphi", "4C")
plotCumFPKMBinArea(tile_width = 10000, all_coverage_counts_df_merged = MT2_all_coverage_counts_df_merged, which_element == "MT2_Mm")

# ORR1A0 plot
ORR1A0_all_coverage_counts_list_merged <- lapply(X = 1:length(all_stages_counts), FUN = function(X) ORR1A0_all_coverage_counts[[X]][, 1:2])
ORR1A0_all_coverage_counts_df_merged <- Reduce(function(...) merge(..., by = 'pos', all.x = TRUE), ORR1A0_all_coverage_counts_list_merged)
colnames(ORR1A0_all_coverage_counts_df_merged) <- c("pos", "GV", "1C", "2C", "2C_aphi", "4C")
plotCumFPKMBinArea(tile_width = 10000, all_coverage_counts_df_merged = ORR1A0_all_coverage_counts_df_merged, which_element == "ORR1A0")


