# seqinr
for(LTR_class_sample in unique(LTR_data$LTR_class)){
  
  aligned <- read.alignment(file = paste0(getwd(), "/", LTR_name, "_random_200", "_aligned.fasta"), 
                            format = "fasta", 
                            forceToLower = F)
  consensus(aligned, method = "majority", type = "DNA")
  
}

# seqlogo
LTR_name <- "MTA"
aligned_seq <- 
  readDNAStringSet(filepath = paste0(getwd(), "/", LTR_name, "_random_200", "_aligned.fasta"), format = 'fasta') %>% 
  consensusMatrix(as.prob = TRUE, baseOnly = T) %>% 
  .[1:4, ]

# msa
library(msa)
aligned_seq_msa <- msa(readDNAStringSet(filepath = paste0(getwd(), "/", LTR_name, "_random_200", "_aligned.fasta"), format = 'fasta'), "ClustalOmega")
msaPrettyPrint(aligned_seq_msa, 
#                y = c(10, 11), 
               subset = c(1:10), 
               showNames = "none", 
               showConsensus = "none",
               showLogo = "bottom", 
               showNumbering = "none",
               shadingMode = "similar",
               showLegend = FALSE, 
               askForOverwrite = FALSE, 
               file = paste0(getwd(), "/", LTR_name, "_random_200_aligned_msaLogo.pdf"), 
               paperWidth = 50,
               paperHeight = 10,
               furtherCode = "\\clearlogocolors \\logocolor{G}{Yellow} \\logocolor{A}{Green} \\logocolor{TU}{Red} \\logocolor{C}{Blue}")

# RWebLogo
library(RWebLogo)
weblogo(file.in = paste0(getwd(), "/", LTR_name, "_random_200", "_aligned.fasta"),
        open = F, 
        file.out = paste0(getwd(), "/", LTR_name, "_random_200", "_aligned_weblogo.pdf"),
        format = "pdf", 
        color.scheme = "classic")
