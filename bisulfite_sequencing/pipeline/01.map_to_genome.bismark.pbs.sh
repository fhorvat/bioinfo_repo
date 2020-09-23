#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.bismark
#PBS -l select=ncpus=20:mem=60g
#PBS -j oe
#PBS -J 0-%N_SAMPLES
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads 
THREADS=4

# input
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# mapping parameters
SCRIPT_PAR="--genome_folder ${BISMARK_INDEX} \
--parallel ${THREADS} \
--unmapped \
--output_dir . \
--temp_dir ."

# ----------------Commands------------------- #
# mapping
if [ ${SINGLE_END} == TRUE ]
then
   bismark ${SCRIPT_PAR} --genome_folder ${BISMARK_INDEX} --single_end ${FILE} 
else
   bismark ${SCRIPT_PAR} --genome_folder ${BISMARK_INDEX} -1 ${FILE}_1.txt.gz -2 ${FILE}_2.txt.gz
fi

# rename
mv ${BASE}.txt.gz_bismark_bt2.bam ${BASE}_bismark_bt2.bam
mv ${BASE}.txt.gz_bismark_bt2_SE_report.txt ${BASE}_bismark_bt2_SE_report.txt
mv ${BASE}.txt.gz_unmapped_reads.fq.gz ${BASE}.unmapped.txt.gz

# sort bam
samtools sort -@ ${THREADS} -o ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2.bam

# rename over original file
[ -f "${BASE}_bismark_bt2.sorted.bam" ] && mv ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2.bam

# bam index
samtools index -@ ${THREADS} ${BASE}_bismark_bt2.bam
