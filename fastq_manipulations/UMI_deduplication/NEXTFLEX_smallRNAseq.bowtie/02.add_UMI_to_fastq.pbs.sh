#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.add_UMI_to_fastq
#PBS -J 0-7
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=..
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.txt.gz" -and -name "*atrim*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.atrim.txt.gz}

# ----------------Commands------------------- #
# ungzip fastq
unpigz -p $THREADS -c $FILE > ${BASE}.fastq

# add the UMI to the fastq file identifier line
awk '{getline p<f} (NR%4==1){$1=$1" "$2;$2=p}1' OFS= f=${BASE}.05.umi_final.txt ${BASE}.fastq > ${BASE}.umi.fastq

# remove reads from fastq with Ns in the UMI:
sed -e '/_N\|_.*N/,+3d' ${BASE}.umi.fastq > ${BASE}.n_removed.fastq

# remove random 4 base pair seqs that make up the UMI from the fastq read sequence line:
cutadapt -u 4 -o ${BASE}.trim_1.fastq ${BASE}.n_removed.fastq
cutadapt -m 18 -u  -4 -o ${BASE}.trim_2.fastq ${BASE}.trim_1.fastq

# remove space form the identifier of the fastq
sed 's/ /-/' ${BASE}.trim_2.fastq > ${BASE}.no_space.fastq

# remove intermediate fastq files
[ -f "${BASE}.no_space.fastq" ] && rm ${BASE}.fastq ${BASE}.umi.fastq ${BASE}.n_removed.fastq ${BASE}.trim_1.fastq ${BASE}.trim_2.fastq 
