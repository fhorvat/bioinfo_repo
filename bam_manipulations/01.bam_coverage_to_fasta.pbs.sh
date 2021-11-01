#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.test_script
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# files
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam" -and -name "*test*" ))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

IN_FASTA=($(find $INPUT_DIR -maxdepth 1 -name "*.fa"))
FILE_FASTA=${IN_FASTA[0]}
BASE_FASTA=${FILE_FASTA#${INPUT_DIR}/}
BASE_FASTA=${BASE_FASTA%.fa}

# ----------------Commands------------------- #
# bam to bedGraph
bedtools genomecov -ibam ${FILE} -bg -split > ${BASE}.bedGraph

# filter by coverage depth
awk '{ if ($4 > 20) { print } }' ${BASE}.bedGraph > ${BASE}.filtered.bedGraph

# merge all overlapping intervals
bedtools merge -i ${BASE}.filtered.bedGraph > ${BASE}.merged.bed

# extract fasta from bed
bedtools getfasta -fi ${FILE_FASTA} -bed ${BASE}.merged.bed > ${BASE}.${BASE_FASTA}.fasta
