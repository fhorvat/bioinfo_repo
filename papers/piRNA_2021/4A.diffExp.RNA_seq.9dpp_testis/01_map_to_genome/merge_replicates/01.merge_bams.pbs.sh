#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_bams
#PBS -l select=ncpus=6:mem=30g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
THREADS=6

STAR_INDEX=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/STAR_index.2.7/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 \( -name "*.bam" -not -name "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE.bam" \)))
STAGES=(${IN_SEQ[@]%_r[1-9]*.[S,P]E.*bam})
STAGES=(${STAGES[@]%_[F,S,M]*})
STAGES=(${STAGES[@]#${INPUT_DIR}/})
STAGES=($(printf "%s\n" ${STAGES[@]} | sort -u))

STAGE=${STAGES[1]}

FILES=($(find $INPUT_DIR -maxdepth 1 \( -name "$STAGE*r*.*E.*bam" -not -name "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE.bam" \)))
FILES=$(printf "%s " ${FILES[@]})

# just in case
echo $FILES

# ----------------Commands------------------- #
# merge
sambamba merge ${STAGE}.bam $FILES -t $THREADS
