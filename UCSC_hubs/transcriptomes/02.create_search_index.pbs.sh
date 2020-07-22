#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.create_search_index
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -name "*.bb" -and \( -name "ensembl*" -or -name "GCF*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bb}

INPUT_DIR=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation
IN_INDEX=($(find ${INPUT_DIR} -name "$BASE*" -and -name "*geneID_to_geneSymbol.txt"))

# ----------------Commands------------------- #
# create index
ixIxx ${IN_INDEX} ${BASE}.ix ${BASE}.ixx

# set permissions
chmod 744 ${BASE}.ix ${BASE}.ixx
