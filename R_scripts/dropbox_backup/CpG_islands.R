CpG.islands.data <- read.table("141107_CpG_islands", header = F, sep = "\t")
colnames(CpG.islands.data) <- c("bin",  "chrom",	"chromStart",	"chromEnd",	
                                "name",	"length",	"cpgNum",	"gcNum",	
                                "perCpg",	"perGc",	"obsExp")
