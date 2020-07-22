zcat GCF_000349665.1_MesAur1.0_rna.fna.gz | grep ">" | sed 's/ .*(/ /g' | sed 's/).*//g' | sed 's/>//g' > GCF_000349665.1_MesAur1.0_rna.geneID_to_geneSymbol.txt
