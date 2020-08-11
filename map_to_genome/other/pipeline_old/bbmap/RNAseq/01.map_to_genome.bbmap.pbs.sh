#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.bbmap
#PBS -l select=ncpus=12:mem=60g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=60G

STAR_INDEX=%GENOME_DIR/STAR_index/%SJDB_OVERHANG
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=../../Raw/Links
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.txt.gz"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# mapping
bbmap.sh -Xmx$MEMORY threads=$THREADS \
in=${FILE}_1.txt.gz \
in2=${FILE}_2.txt.gz \
outm=${BASE}.sam outu=${BASE}.unmapped.fq.gz statsfile=${BASE}.log.out \
perfectmode=t ambiguous=all pairedonly=t maxsites=20000 \
nhtag=t amtag=t scoretag=t unpigz=t trimreaddescriptions=t

# sort sam and transform to bam, index
samtools view -bShu ${BASE}.sam | samtools sort -@ $THREADS - -o ${BASE}.bam
samtools index ${BASE}.bam

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph

# permissions
chmod 744 ${BASE}.bam ${BASE}.bam.bai ${BASE}.bw

