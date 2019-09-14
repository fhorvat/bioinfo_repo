#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.map_to_rDNA_45S_mouse.STAR
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4
STAR_INDEX=%GENOME_DIR/STAR_index/rDNA_45S_BK000964.3

INPUT_DIR=.
IN_SEQ=(`ls $INPUT_DIR/*.genome.Unmapped.out.mate1`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.genome.Unmapped.out.mate1}.rDNA_45S

STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM  20000000000 \
--outFileNamePrefix ${BASE}. \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 99999 \
--outFilterMultimapScoreRange 0 \
--outFilterMismatchNoverLmax 0.2"

# ----------------Commands------------------- #
# mapping
STAR --readFilesIn ${FILE} $STAR_PAR
