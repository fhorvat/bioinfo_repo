#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.add_UMI_to_fastq
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.filtered.fastq" -not -name "*umi*" -not -name "*n_removed*" -not -name "*trim*" -not -name "*no_space*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.filtered.fastq}

# ----------------Commands------------------- #
# add the UMI to the fastq file identifier line
awk '{getline p<f} (NR%4==1){$1=$1""$2;$2=p}1' OFS= f=${BASE}.05.umi_final.txt ${FILE} > ${BASE}.umi.fastq

# remove reads from fastq with Ns in the UMI:
sed -e '/_N\|_.*N/,+3d' ${BASE}.umi.fastq > ${BASE}.n_removed.fastq

# remove random 4 base pair seqs that make up the UMI from the fastq read sequence line:
cutadapt -u 4 -o ${BASE}.trim_1.fastq -j $THREADS ${BASE}.n_removed.fastq
cutadapt -m 18 -u  -4 -o ${BASE}.trim_2.fastq -j $THREADS ${BASE}.trim_1.fastq

# remove space from the identifier of the fastq
#sed 's/ /-/' ${BASE}.trim_2.fastq > ${BASE}.no_space.fastq

# remove intermediate fastq files
[ -f "${BASE}.trim_2.fastq" ] && rm ${BASE}.umi.fastq ${BASE}.n_removed.fastq ${BASE}.trim_1.fastq

# rename
mv ${BASE}.trim_2.fastq ${BASE}.fastq

# gzip
pigz -p $THREADS ${BASE}.fastq

# rename
mv ${BASE}.fastq.gz ${BASE}.txt.gz
