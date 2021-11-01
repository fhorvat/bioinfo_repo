#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bam_to_stranded_bw
#PBS -l select=ncpus=1:mem=30g
#PBS -J 0-11
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
THREADS=1

STAR_INDEX=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/STAR_index/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

SCRIPT=/common/WORK/fhorvat/programi/stranded_coverage/stranded-coverage/strand_cov

# ----------------Commands------------------- #
# stranded coverage from bam
$SCRIPT -o $BASE $FILE

# convert to bigwig
wigToBigWig ${BASE}.plus.wig $CHR_LENGTH ${BASE}.plus.bw
wigToBigWig ${BASE}.minus.wig $CHR_LENGTH ${BASE}.minus.bw

# remove wigs
[ -f "${BASE}.plus.bw" ] && rm ${BASE}.plus.wig
[ -f "${BASE}.minus.bw" ] && rm ${BASE}.minus.wig
