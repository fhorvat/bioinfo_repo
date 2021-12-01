library(rentrez)
r_search <- entrez_search(db = "sra",
                          term = "SRP020490[ACCN] AND 4-cell[WORD]",
                          retmax = 20)

entrez_sra_summary <- entrez_summary(db = "sra", r_search$ids)
entrez_sra_summary <- extract_from_esummary(entrez_sra_summary, "runs")
entrez_sra_summary <- unlist(strsplit(entrez_sra_summary, split = " "))
entrez_sra_summary <- grep("^acc=", entrez_sra_summary, value = T)
entrez_sra_summary <- gsub("^acc=\\\"|\\\"$", "", entrez_sra_summary)
