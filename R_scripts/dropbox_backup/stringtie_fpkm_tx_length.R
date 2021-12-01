library("rtracklayer")
library("dplyr")
library("reshape")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MousePolyA/2cell/Hamazaki_2015")

filenames <- list.files(path = "mapping/STAR/StringTie", pattern = "\\.gtf$", full.names = T)
sample_names <- list.files(path = "mapping/STAR/StringTie", pattern = "\\.gtf$", full.names = F)
sample_names <- gsub("_.*", "", sample_names)

# list of .gtf tables
gtf_tables <- lapply(filenames, import)

# transcript width (sum of lengths of exons which belong to each transcript)
tx_width <- data.frame(type = mcols(gtf_tables[[1]])$type, 
                       width = width(gtf_tables[[1]]), 
                       gene_id = mcols(gtf_tables[[1]])$gene_id,
                       transcript_id = mcols(gtf_tables[[1]])$transcript_id, 
                       gene_name = mcols(gtf_tables[[1]])$ref_gene_name)
tx_width <- group_by(tx_width, transcript_id, type == "exon")

tx_width_sum <- summarise(tx_width, sum(width))
colnames(tx_width_sum) <- c("transcript_id", "type", "width_sum")
tx_width_sum <- filter(tx_width_sum, type == "TRUE")

tx_width <- full_join(tx_width, tx_width_sum, by = "transcript_id")
tx_width <- filter(tx_width, type.x == "transcript")
tx_width <- data.frame(tx_width)
tx_width <- tx_width[, c("transcript_id", "width_sum", "gene_id", "gene_name")]


# list of FPMK for each sample in one data.frame
FPMK_df <- lapply(X = 1:length(gtf_tables),
              function(X){
                data.frame(transcript_id = mcols(gtf_tables[[X]])$transcript_id, 
                           FPKM = as.numeric(mcols(gtf_tables[[X]])$FPKM),
                           type = mcols(gtf_tables[[X]])$type)
              })
FPMK_df <- lapply(FPMK_df, filter, type == "transcript")
FPMK_df <- lapply(X = FPMK_df, function(X) X[, 1:2])
names(FPMK_df) <- sample_names

lapply(X = 1:length(FPMK_df), 
       function(X){
         colnames(FPMK_df[[X]]) <<- c("transcript_id", 
                                  paste0("FPKM_", names(FPMK_df)[X]))
       })
FPMK_df <- Reduce(function(...) merge(..., by = "transcript_id", all = T), FPMK_df)

# merge FPMKs for each sample with transcript width
FPMK_tx_width <- merge(tx_width, FPMK_df, by = "transcript_id", all = T)

# remove transcript with FPKM = 0 in all samples
# FPMK_tx_width <- FPMK_tx_width[rowSums(FPMK_tx_width[ ,names(FPMK_tx_width)[grep("FPKM", names(FPMK_tx_width))]]) > 0, ]

write.csv(FPMK_tx_width, "FPMK_tx_width.csv", row.names = F)
