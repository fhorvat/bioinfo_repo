#/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.bowtie2
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-20
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
BOWTIE2_INDEX_NAME=myco

INPUT_DIR=../../Raw/Links
IN_SEQ=(`ls ${INPUT_DIR}/*.txt.gz`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}
 
# print number of samples
echo "Number of samples:" ${#IN_SEQ[@]}

# set script parameters
BOWTIE2_PAR="--sensitive \
--local \
--no-unal \
--threads $THREADS"

# ----------------Commands------------------- #
# mapping
bowtie2 $BOWTIE2_PAR -x $BOWTIE2_INDEX_NAME -U $FILE -S ${BASE}.sam 2> ${BASE}.log

# sam to bam
samtools view ${BASE}.sam -@ $THREADS -Sb > ${BASE}.bam

# sort bam
samtools sort -@ $THREADS -o ${BASE}_sorted.bam ${BASE}.bam
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.bam
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.sam
mv ${BASE}_sorted.bam ${BASE}.bam

# index
samtools index ${BASE}.bam
