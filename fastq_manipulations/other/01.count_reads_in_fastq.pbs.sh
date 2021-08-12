#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
<<<<<<< HEAD
#PBS -N pbs.01.countFastq
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-3
=======
#PBS -N pbs.countFastq
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-15
>>>>>>> 15478517ec9c5615c16fe81bcc2371bd062724d9
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
<<<<<<< HEAD
IN_SEQ=($(find ${INPUT_DIR} -name "*.fastq.gz" -or -name "*.txt.gz"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fastq.gz}
=======
IN_SEQ=($INPUT_DIR/*txt.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
>>>>>>> 15478517ec9c5615c16fe81bcc2371bd062724d9
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# counts reads in fastq
<<<<<<< HEAD
echo ${BASE} `zcat ${FILE} | awk '{s++}END{print s/4}'` >> stats.txt
=======
echo ${BASE} `zcat ${BASE}.txt.gz | awk '{s++}END{print s/4}'` >> stats.txt
>>>>>>> 15478517ec9c5615c16fe81bcc2371bd062724d9
# sort -t1 stats.txt > stats_sort.txt
