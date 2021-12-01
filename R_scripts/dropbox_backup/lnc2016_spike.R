# lnc2016
sample_table_lnc2016_full <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/RNAseq_2016_11_23_sampleTable.csv", header = T) %>%
  dplyr::select(ID = sample, Time.Course = Sample, Treatment.Control = Type) %>%
  mutate(Treatment.Control = gsub(" B6", "", Treatment.Control), 
         ID = gsub("_16.*", "", ID), 
         name = paste(gsub(".*_[1-8]_", "", ID), Treatment.Control, sep = "_"), 
         experiment = "lnc2016") %>% 
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
  mutate(name = paste(ID, Time.Course, Treatment.Control, sep = "_")) %>%
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

sample_table_lnc2016_full[, c("ID", "Treatment.Control", "lib_size", "lib_size_spike", "ratio")]
sample_table_CNOT6L[, c("ID", "Time.Course", "Treatment.Control", "lib_size", "lib_size_spike", "ratio")]

write.table(file = "CNOT6L_lib_spike.txt", x = sample_table_CNOT6L[, c("ID", "Time.Course", "Treatment.Control", "lib_size", "lib_size_spike")], sep = "\t")
