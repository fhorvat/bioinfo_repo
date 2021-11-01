#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.add_header
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set track name
NAME=PRDM9_ChIP

# input files
INPUT_DIR=.
IN_SEQ=($(find $INPUT_PATH -maxdepth 1 -name "*narrowPeak" -not -name "*UCSC*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.narrowPeak}

# ----------------Commands------------------- #
# create new .narrowPeak file with header for UCSC genome browser
echo "track type=narrowPeak name="$NAME"" > ${BASE}.UCSC.narrowPeak
cat $FILE >> ${BASE}.UCSC.narrowPeak

# set permission
chmod 744 ${BASE}.UCSC.narrowPeak
