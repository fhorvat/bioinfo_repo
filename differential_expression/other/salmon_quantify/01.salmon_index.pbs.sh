#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.salmon_index
#PBS -l select=ncpus=4:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4

INPUT_DIR=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1
IN_SEQ=(${INPUT_DIR}/ensembl.93.MesAur1.0.20180920.cdna.all.fa.gz)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE/.fa.gz/.salmon_index}

# ----------------Commands------------------- #
# make salmon index from ENSEMBL cDNA fasta
salmon index -t ${FILE} -i ${BASE} -p $THREADS
