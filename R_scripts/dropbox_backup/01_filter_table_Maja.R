setwd("C:/Users/fhorvat/Dropbox/Praksa bioinfo/DESeq_analysis_Svoboda/Prague/MT_elements/MT2_solo_ORR1A0_coverage")
load("mt2_mm_orr1a0.Robj")

mt2_mm_orr1a0_df <- as.data.frame(mt2_mm_orr1a0)
mt2_mm_orr1a0_df <- mt2_mm_orr1a0_df[, c(1, 2, 3, 5, 6)]

mt2_mm_orr1a0_df_ORR <- mt2_mm_orr1a0_df[mt2_mm_orr1a0_df$repName == "ORR1A0", ]
mt2_mm_orr1a0_df_MT2 <- mt2_mm_orr1a0_df[mt2_mm_orr1a0_df$repName == "MT2_Mm", ]

write.table(mt2_mm_orr1a0_df_ORR, 
            "ORR1A0_solo_LTRs_Maja_20160803.txt", 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

write.table(mt2_mm_orr1a0_df_MT2, 
            "MT2_solo_LTRs_Maja_20160803.txt", 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)
