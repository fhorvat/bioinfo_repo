library("dplyr")

rptmsk_org <- read.delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20160823.fa.out.txt.gz", stringsAsFactors = F, header = T)

rptmsk_viz <- read.delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz", stringsAsFactors = F, header = T)
rptmsk_viz <- rptmsk_viz[, c("genoName", "genoStart", "genoEnd", "strand", "repName", "repClass", "id")]
colnames(rptmsk_viz) <- c("seqnames", "start", "end", "strand", "element_name", "element_class", "VIZ_ID")

# filter for MT, ORR and MLT
rptmsk_org_filtered <- rptmsk_org[grep("SINE", rptmsk_org$element_class), ]
rptmsk_viz_filtered <- rptmsk_viz[grep("SINE", rptmsk_viz$element_class), ]

# number of fragments repeat masker
rptmsk_org_B2 <- rptmsk_org_filtered[grep("B2", rptmsk_org_filtered$element_class), ]

rptmsk_viz_B2 <- rptmsk_viz_filtered[grep("B2|B3", rptmsk_viz_filtered$element_name), ]
rptmsk_viz_B2_VIZ_ID_unique <-
  rptmsk_viz_B2 %>%
  group_by(element_name) %>%
  summarise(element_name_number = length(unique(VIZ_ID)))

