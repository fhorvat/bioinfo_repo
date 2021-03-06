#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03b.map_to_rDNA_45S_mouse.STAR
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4
GENOME_DIR=%GENOME_DIR/STAR_index/rDNA_45S_BK000964.3

INPUT_DIR=.
IN_SEQ=(`ls $INPUT_DIR/*.genome.merged.Unmapped.out.mate1`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.genome.merged.Unmapped.out.mate1*}.rDNA_45S

STAR_PAR="--genomeDir $GENOME_DIR \
--runThreadN $THREADS \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM  20000000000 \
--outFileNamePrefix ${BASE}. \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 99999 \
--outFilterMismatchNoverLmax 0.2"

# ----------------Commands------------------- #
# mapping
STAR --readFilesIn ${FILE} $STAR_PAR
