all_bz - "((rn6 (GraSur mm10)) (AcoCah MerUng))" >&all_bz.log
sed -i 's/Y=9000 H=0/M=40 K=4500 L=4500 Y=15000 T=2 O=600 E=55/' all_bz.log

sh all_bz.log

tba "((rn6 (GraSur mm10)) (AcoCah MerUng))" *.*.maf tba.maf >&tba.log

# create bigBed tracks
GENOMES=(AcoCah MerUng GraSur)
for GENOME in ${GENOMES[@]} 
do
	maf_project tba.maf $GENOME > $GENOME.maf
	mafToBigMaf $GENOME $GENOME.maf stdout | sort -k1,1 -k2,2n > $GENOME.bigMaf
	bedToBigBed -type=bed3+1 -as=/common/WORK/fhorvat/programi/UCSC/bigMaf.as -tab $GENOME.bigMaf ../../../chrom_sizes/$GENOME*chrom.sizes $GENOME.bb
done
