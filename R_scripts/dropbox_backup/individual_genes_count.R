count_df <- read.csv(file = "/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Analysis/counts_knownGene_mm10_lncNov2016_CNOT6L.csv", header = T, stringsAsFactors = F)
colnames(count_df)[1] <- "entrezID"
count_df <- count_df[, grepl("GV_WT|GV_KO|BC|entrezID", colnames(count_df))]
count_df$gene_symbol <- mapIds(org.Mm.eg.db,
                               keys = as.character(count_df$entrezID),
                               column ="SYMBOL",
                               keytype = "ENTREZID",
                               multiVals = "first")
count_df$gene_name <- mapIds(org.Mm.eg.db,
                               keys = as.character(count_df$entrezID),
                               column ="GENENAME",
                               keytype = "ENTREZID",
                               multiVals = "first")

# lnc2016
sample_table_lnc2016_full <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/RNAseq_2016_11_23_sampleTable.csv", header = T) %>%
  dplyr::select(ID = sample, Time.Course = Sample, Treatment.Control = Type) %>%
  mutate(Treatment.Control = gsub(" B6", "", Treatment.Control), 
         ID = gsub("_16.*", "", ID), 
         name = paste(gsub(".*_[1-8]_", "", ID), Treatment.Control, sep = "_")) %>% 
  left_join(data.frame(sample_paths = list.files("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                                 pattern = "*.bam$",
                                                 recursive = T,
                                                 full.names = T), 
                       logs_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                               pattern = "*Log.final.out", 
                                               recursive = T, 
                                               full.names = T), 
                       spike_log_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/spike_mapping/mapped",
                                                    pattern = "*Log.final.out",
                                                    recursive = T,
                                                    full.names = T),
                       stringsAsFactors = F) %>%
              mutate(ID = gsub("^/.*/|_16.*", "", sample_paths)) %>%
              rowwise() %>%
              mutate(lib_size = (as.integer(read.delim(logs_paths, header = F, stringsAsFactors = F)[8, 2])), 
                     lib_size_spike = (as.integer(read.delim(spike_log_paths, header = F, stringsAsFactors = F)[8, 2]))), 
            by = "ID")

# CNOT6L
sample_table_CNOT6L <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_sample_list_11919R_2015_10_29.csv", header = T) %>%
  dplyr::select(ID, Time.Course, Treatment.Control) %>%
  mutate(name = paste(paste0("X", ID), Time.Course, Treatment.Control, sep = "_")) %>%
  left_join(data.frame(sample_paths = list.files("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                                 pattern = "*.bam$",
                                                 recursive = T, 
                                                 full.names = T), 
                       logs_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                               pattern = "*Log.final.out", 
                                               recursive = T, 
                                               full.names = T), 
                       spike_log_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/spike_mapping/mapped",
                                                    pattern = "*Log.final.out",
                                                    recursive = T,
                                                    full.names = T), 
                       stringsAsFactors = F) %>%
              mutate(ID = gsub("^/.*/|_.*", "", sample_paths)) %>%
              rowwise() %>%
              mutate(lib_size = (as.integer(read.delim(logs_paths, header = F, stringsAsFactors = F)[8, 2])), 
                     lib_size_spike = (as.integer(read.delim(spike_log_paths, header = F, stringsAsFactors = F)[8, 2]))), 
            by = "ID")

table_full <- rbind(sample_table_lnc2016_full, sample_table_CNOT6L)

df <- as.data.frame(t(count_df[grepl("Zp3$", count_df$gene_symbol) | grepl("Mos$", count_df$gene_symbol), ]))
colnames(df) <- c("Mos", "Zp3")
df <- df[!grepl("entrezID|gene_symbol|gene_name", rownames(df)), ]
df$name <- rownames(df)

df <- 
  left_join(df, table_full, by = "name") %>%
  select(name, Mos, Zp3, lib_size, lib_size_spike )

write.csv(x = df, file = "Zp3_Mos_GV_counts.csv")
