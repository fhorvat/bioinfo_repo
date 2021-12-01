library(tabulizer)

# Extract the table
out <- extract_tables("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4434936/bin/supp_112.105312_105312SupData.pdf")

pdf_df <- 
  lapply(out, as.data.frame) %>% 
  dplyr::rbind_all(.)