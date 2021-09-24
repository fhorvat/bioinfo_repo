#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.UMI_list
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.atrim.txt.gz" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.atrim.txt.gz}

# ----------------Commands------------------- #
# reading every 4th line starting with line 2, get first 4 characters of sequence
zcat ${FILE} | awk 'NR%4==2' | cut -d' ' -f2 | cut -c1-4 > ${BASE}.01.umi_first.txt

# reading every 4th line starting with line 2, get last 4 characters of sequence
zcat ${FILE} | awk 'NR%4==2' | sed 's/^.*\(.\{4\}\)/\1/' > ${BASE}.02.umi_last.txt

# pasting first UMI 4 nuc. with last UMI 4 nuc.
paste -d'\0' ${BASE}.01.umi_first.txt ${BASE}.02.umi_last.txt > ${BASE}.03.umi_together.txt

# quadruple UMIs
awk '{for(i=0;i<4;i++)print}' ${BASE}.03.umi_together.txt > ${BASE}.04.umi_quad.txt

# add an "_" to the front of every UMI line
awk '$0="_"$0' ${BASE}.04.umi_quad.txt > ${BASE}.05.umi_final.txt

# remove intermediate files
[ -f "${BASE}.05.umi_final.txt" ] && rm ${BASE}.01.umi_first.txt ${BASE}.02.umi_last.txt ${BASE}.03.umi_together.txt ${BASE}.04.umi_quad.txt
