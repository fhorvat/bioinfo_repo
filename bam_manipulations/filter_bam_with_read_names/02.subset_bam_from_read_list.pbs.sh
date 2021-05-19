#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.filter_bam
#PBS -l select=ncpus=6:mem=50g
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=6
MEMORY=50g

# files
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.deduplicated.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam" -and -name "*18to32nt*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

INPUT_DIR=../bam_subset
NAMES_FILE=${INPUT_DIR}/${BASE}.read_list.txt

# ----------------Commands------------------- #
### extract all alignments of reads provided in a list of read names
# awk command first parses the reads.txt file, creates an array of QNAMEs,
# then parses the alignments and print if the QNAME is present in the array (or if the line belongs to the header).
samtools view -h ${FILE} \
| awk 'FNR==NR {reads[$1];next} /^@/||($1 in reads)' ${NAMES_FILE} - \
| samtools view -b - > ${BASE}.bam
