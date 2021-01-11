#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.substitution_rate
#PBS -l select=ncpus=1:mem=1g
#PBS -J 0-625
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# run parameters
THREADS=1

# script path
SCRIPT=02b.substitution_rate.script.R

# input for script
INPUT_DIR=./fasta_files
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.single.fasta"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.200730.single.fasta}
BASE=${BASE#LTRs.random_200_per_repName.}

# input for manual script
echo -e '\n'\
fasta_file=\'$FILE\''\n'\
fasta_name=\'$BASE\''\n'\
threads=as.numeric\(\'$THREADS\'\)'\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--fasta_file $FILE \
--fasta_name $BASE \
--threads $THREADS
