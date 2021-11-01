#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N opossum
#PBS -J 0-4
#PBS -l select=ncpus=1:mem=20g
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
OPOSSUM_PATH=/common/WORK/fhorvat/programi/Opossum/Opossum.py

INPUT_DIR=./STAR_2pass
IN_SEQ=(`find $INPUT_DIR -name "WT*.bam" -o -name "Lnc5*.bam"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_Aligned.sortedByCoord.out.bam}

# ----------------Commands------------------- #
# process for variant calling
python $OPOSSUM_PATH --BamFile $FILE --OutFile ${BASE}.bam --SoftClipsExist True --ProperlyPaired False 2> ${BASE}.log
