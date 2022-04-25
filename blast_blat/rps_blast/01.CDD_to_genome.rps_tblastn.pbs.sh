#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.CDD_to_genome.rps_blast
#PBS -l select=ncpus=18:mem=100g
#PBS -j oe
 
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=18

# database file
INPUT_DIR=/common/DB/genome_reference/other/CDD/Cog_LE
IN_DB=${INPUT_DIR}/Cog

# genome file
INPUT_DIR=/common/DB/genome_reference/Mollusca/Arion_vulgaris.Schrodl
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "AriVul.fix.fa"))
FILE_GENOME=${IN_GENOME[0]}
BASE=${FILE_GENOME#${INPUT_DIR}/}
BASE=${BASE%.fa}

# ----------------Commands------------------- #
# blast
rpstblastn -query ${FILE_GENOME} -db ${IN_DB} -out ${BASE}.cdd.txt -num_threads ${THREADS} -outfmt 6

#-gapopen 10 \
#-gapextend 20 \
#-penalty -1
