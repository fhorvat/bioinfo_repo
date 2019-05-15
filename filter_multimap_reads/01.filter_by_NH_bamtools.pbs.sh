#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N bamtools_filter
#PBS -l select=ncpus=2:mem=10g
#PBS -J 0-10
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter
IN_SEQ=($INPUT_DIR/*.bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_Aligned.sortedByCoord.out.bam}

# ----------------Commands------------------- #
# bamtools filter by NH tag
bamtools filter -in $FILE -out ${BASE}_multimap.bam -script multimap_filter.json
