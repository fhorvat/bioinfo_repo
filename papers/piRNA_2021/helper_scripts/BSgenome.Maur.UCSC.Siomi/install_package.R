# instructions from here: https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf

library(BSgenome)
forgeBSgenomeDataPkg("/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/BSgenome.Maur.UCSC.Siomi/BSgenome.Maur.UCSC.Siomi.seed")

# after done run in command line:
# R CMD build BSgenome.Maur.UCSC.Siomi/
# R CMD check BSgenome.Maur.UCSC.Siomi_1.0.0.tar.gz
# R CMD INSTALL BSgenome.Maur.UCSC.Siomi_1.0.0.tar.gz
