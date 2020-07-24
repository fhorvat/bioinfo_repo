#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.bedGraph_to_bigWig
#PBS -l select=ncpus=1:mem=30g
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index.2.7/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bdg"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bdg}

# ----------------Commands------------------- #
# bedGraph to bigWig
wigToBigWig ${BASE}.bdg ${CHR_LENGTH} ${BASE}.bw
#[ -f "${BASE}.bw" ] && rm ${BASE}.bdg

# set permissions
chmod 744 ${BASE}.bw
