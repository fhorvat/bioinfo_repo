library(SplicingGraphs)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Gviz)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

isActiveSeq(txdb)[1:25]

isActiveSeq(txdb)[match("chr14", names(isActiveSeq(txdb)))] <- TRUE
names(which(isActiveSeq(txdb)))

sg <- SplicingGraphs(txdb)
str(sg)
sg[["3183"]]
mcols(sg[["3183"]])
mcols(sg[["3183"]])$txpath

plotTranscripts(sg[["3183"]])

ex_by_tx <- unlist(sg)
ex_by_tx[names(ex_by_tx) %in% c("10001", "100129075")]

sg[strand(sg) == "-"]
sg[1:20]

tail(sg)

sgedges(sg["3183"])
edges_by_gene <- sgedgesByGene(sg)

plot(sg["3183"])
plot(sgraph(sg["3183"], tx_id.as.edge.label=TRUE))


AScode_list <- lapply(seq_along(sg), function(i) bubbles(sg[i])$AScode)
names(AScode_list) <- names(sg)
AScode_table <- table(unlist(AScode_list))
AScode_table <- sort(AScode_table, decreasing=TRUE)

AScode_summary <- data.frame(AScode=names(AScode_table),
                             NbOfEvents=as.vector(AScode_table),
                             Desciption=ASCODE2DESC[names(AScode_table)])


