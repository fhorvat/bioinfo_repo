#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=5g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.fix_CDS_in_Liftoff_output.gt
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=5G

# query .gtf
INPUT_DIR=.
IN_GTF=($(find $INPUT_DIR -maxdepth 1 \( -name "*.liftoff.gff" -not -name "*fixed.gff" \)))
FILE_GTF=${IN_GTF[0]}
BASE_GTF=${FILE_GTF#${INPUT_DIR}/}
BASE_GTF=${BASE_GTF%.gff}

# script
SCRIPT=/common/WORK/fhorvat/programi/genometools/gt-1.5.10-Linux_x86_64-64bit-barebone/bin/gt

# ----------------Commands------------------- #
# liftoff in versions up to 1.61 doesn't assign correct phase to CDS
# this process is used to correct this:
$SCRIPT gff3 -sort -tidy -retainids ${FILE} > ${BASE}.fixed.gff
