#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_bam
#PBS -l select=ncpus=12:mem=50g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
STAR_INDEX=/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6/STAR_index.2.7/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=.
IN_SEQ=$(find $INPUT_DIR \( -name "*bam" -not -name "*merged.bam" \) | paste -sd" " -)
BASE="s_GoldHam_all.spliced.merged"

# ----------------Commands------------------- #
# merge
sambamba merge -t $THREADS ${BASE}.bam ${IN_SEQ}
