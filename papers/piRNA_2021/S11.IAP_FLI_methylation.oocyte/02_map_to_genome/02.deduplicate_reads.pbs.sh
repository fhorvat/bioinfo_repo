#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.deduplicate_reads
#PBS -l select=ncpus=16:mem=20g
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads
THREADS=4
MEMORY=20G

# input
IN_SEQ=($(find . -maxdepth 1 \( -name "*.bam" -not -name "*deduplicated.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# run 
deduplicate_bismark --outfile ${BASE} --output_dir . ${FILE}

# get coverage
genomeCoverageBed -ibam ${BASE}.deduplicated.bam -bg -split > ${BASE}.deduplicated.raw.bedGraph
wigToBigWig ${BASE}.deduplicated.raw.bedGraph ${CHR_LENGTH} ${BASE}.deduplicated.raw.bw
[ -f "${BASE}.deduplicated.raw.bw" ] && rm ${BASE}.deduplicated.raw.bedGraph
