#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.bismark
#PBS -l select=ncpus=10:mem=20g
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads 
THREADS=4

# input
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[0,1,2,s].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# mapping parameters
SCRIPT_PAR="--genome_folder ${BISMARK_INDEX} \
--non_directional \
--parallel ${THREADS} \
--output_dir . \
--temp_dir ."

# ----------------Commands------------------- #
# mapping
if [ ${SINGLE_END} == TRUE ]
then

   # run bismark
   bismark ${SCRIPT_PAR} --genome_folder ${BISMARK_INDEX} --single_end ${FILE}.fastq

   # sort bam
   samtools sort -@ ${THREADS} -o ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2.bam

   # rename over original file
   [ -f "${BASE}_bismark_bt2.sorted.bam" ] && mv ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2.bam

   # bam index
   samtools index -@ ${THREADS} ${BASE}_bismark_bt2.bam   

else

   # run bismark
   bismark ${SCRIPT_PAR} --genome_folder ${BISMARK_INDEX} -1 ${FILE}_1.txt.gz -2 ${FILE}_2.txt.gz ${FILE}_s.txt.gz

   # rename
   mv ${BASE}_1.txt.gz_bismark_bt2_pe.bam ${BASE}_bismark_bt2_pe.bam
   mv ${BASE}_1.txt.gz_bismark_bt2_PE_report.txt ${BASE}_bismark_bt2_PE_report.txt
   mv ${BASE}_1.txt.gz_unmapped_reads_1.fq.gz ${BASE}_1.unmapped.txt.gz
   mv ${BASE}_2.txt.gz_unmapped_reads_2.fq.gz ${BASE}_2.unmapped.txt.gz

   # sort bam
   samtools sort -@ ${THREADS} -o ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2_pe.bam

   # rename over original file
   [ -f "${BASE}_bismark_bt2.sorted.bam" ] && mv ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2_pe.bam

   # bam index
   samtools index -@ ${THREADS} ${BASE}_bismark_bt2_pe.bam

fi
