#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.bismark
#PBS -l select=ncpus=6:mem=30g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads 
THREADS=6

# input
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[1,2,s].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# mapping parameters
SCRIPT_PAR="-OT ${THREADS} \
-UN \
-CP 0.9 \
-L 100,100 \
-t ${THREADS}"

# ----------------Commands------------------- #
# mapping
if [ ${SINGLE_END} == TRUE ]
then

   # run bsbolt
   bsbolt Align -F1 ${FILE}.txt.gz -DB ${BSBOLT_INDEX} -O ${BASE} ${SCRIPT_PAR}

   # bam index
   #samtools index -@ ${THREADS} ${BASE}.bam

   # bam to scaled bedGraph, bedGraph to bigWig
   #genomeCoverageBed -ibam ${BASE}.bam -bg -split > ${BASE}.bedGraph
   #wigToBigWig ${BASE}.bedGraph ${CHR_LENGTH} ${BASE}.raw.bw
   #[ -f "${BASE}.raw.bw" ] && rm ${BASE}.bedGraph

else
   
   # run bsbolt
   bsbolt Align -F1 ${FILE}_1.txt.gz -F2 ${FILE}_2.txt.gz -DB ${BSBOLT_INDEX} -O ${BASE} ${SCRIPT_PAR}

fi
