#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.CDD_to_proteins.rps_blast
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
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/Mollusca_project/protein_annotation/Arion_vulgaris/protein_sequences_from_gff
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "AriVul.CDS.joined.translated_into_protein.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE=${FILE_GENOME#${INPUT_DIR}/}
BASE=${BASE%.fa}

# ----------------Commands------------------- #
# blast
rpsblast -query ${FILE_GENOME} -db ${IN_DB} -out ${BASE}.cdd.txt -num_threads ${THREADS} -outfmt 6

#-gapopen 10 \
#-gapextend 20 \
#-penalty -1
