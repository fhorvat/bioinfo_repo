#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=5g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.insert_size_bam.Picard
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=5G

IN_DIR=.
IN_SEQ=(${IN_DIR}/*bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.bam}

FASTA_FILE=${IN_DIR}/MesAur1.0.fasta

SCRIPT=/common/WORK/fhorvat/programi/Picard/picard.jar

# ----------------Commands------------------- #
# insert size Picard tools
java -jar $SCRIPT CollectInsertSizeMetrics \
R=${FASTA_FILE} \
I=$FILE \
O=${BASE}.insertSize.txt \
HISTOGRAM_FILE=${BASE}.insertSizeHist.pdf
