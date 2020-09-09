#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.PBAT.trim_galore
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -name "*.txt.gz" -and -name "*PE_*.txt.gz"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[0]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# trim
trim_galore --paired ${FILE}_1.txt.gz ${FILE}_2.txt.gz

# rename
mv ${BASE}_1.txt.gz_trimmed.fq.gz ${BASE}_1.txt.gz
mv ${BASE}_2.txt.gz_trimmed.fq.gz ${BASE}_2.txt.gz
