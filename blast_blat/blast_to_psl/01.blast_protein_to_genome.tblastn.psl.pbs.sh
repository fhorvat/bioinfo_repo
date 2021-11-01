#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.blastp.pbs.sh
#PBS -l select=ncpus=6:mem=30g
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=6

# genome file
INPUT_DIR=/common/DB/genome_reference/Mollusca/Arion_vulgaris.Schrodl
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "AriVul.fix.fa"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fa}
IN_DB=${INPUT_DIR}/blast_database/nucleotide/${BASE_GENOME}

# blastn query input
INPUT_DIR=../nematostella_proteins
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.protein.fa"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

# ----------------Commands------------------- #
# blast
tblastn -query ${FILE} -db ${IN_DB} -out ${BASE}.tblastn.xml -num_threads ${THREADS} -outfmt 5

# blast result to psl
blastXmlToPsl -pslx -tName=Hit_def0 -qName=query-def0 ${BASE}.tblastn.xml ${BASE}.tblastn.psl

# psl score
pslScore ${BASE}.tblastn.psl > ${BASE}.tblastn.score.psl
