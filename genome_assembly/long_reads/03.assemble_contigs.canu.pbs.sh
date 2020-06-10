#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.assemble_contigs.canu
#PBS -l select=ncpus=1:mem=5g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -maxdepth 1 -name "*fasta" -not -name "*assembled*"`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

OUT_DIR=./${BASE}
SCRIPT=/common/WORK/kristian/bin/Canu/canu/Linux-amd64/bin/canu

# ----------------Commands------------------- #
# assemble
$SCRIPT -assemble -p ${BASE} -d ${OUT_DIR} \
genomeSize=2000k \
stopOnReadQuality=false \
minOverlapLength=100 \
minReadLength=100 \
-pacbio-corrected ${FILE}

