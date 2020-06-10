#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N solid2fastq
#PBS -J 0-7
#PBS -l select=ncpus=1:mem=10g
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/Cao_2014_pig/Data/Raw/Links
IN_SEQ=($INPUT_DIR/*.csfasta.gz)
FILE_SEQ=${IN_SEQ[$PBS_ARRAY_INDEX]}
FILE_QUAL=${FILE_SEQ%.csfasta.gz}_QV.qual.gz

BASE=${FILE_SEQ##*$INPUT_DIR/}
BASE=${BASE%_F3.csfasta.gz}

# ----------------Commands------------------- #
# solid2fastq
/common/WORK/fhorvat/programi/bfast/BFAST-bfast.0.6.4a/bin/solid2fastq -o $BASE -z -Z $FILE_SEQ $FILE_QUAL
