#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.deduplex_mapped.UMI_tools
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.mapped.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.mapped.bam}

# ----------------Commands------------------- #
# deduplex using UMItools
umi_tools dedup \
--method=directional --multimapping-detection-method=NH \
-I ${FILE} --output-stats=${BASE}.dedup_stats --log=${BASE}.dedup_log.txt -S ${BASE}.dedup.bam
