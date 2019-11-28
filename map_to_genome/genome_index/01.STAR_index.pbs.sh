#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=10:mem=200g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.STAR_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=10

# set input variables
IN_GENOME=($(find .. -name "*.fa.gz"))
IN_GTF=($(find .. -name "*.gtf.gz"))
SJDB_OVERHANG=249
OUTDIR=./sjdbOverhang_${SJDB_OVERHANG}

# create outdir
mkdir $OUTDIR 2> /dev/null

# ----------------Commands------------------- #
# unzip genome and .gtf
unpigz -p $THREADS ${IN_GENOME}
unpigz -p $THREADS ${IN_GTF}

# generate index
STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $OUTDIR --genomeFastaFiles ${IN_GENOME%.gz} --sjdbGTFfile ${IN_GTF%.gz} --sjdbOverhang $SJDB_OVERHANG --limitGenomeGenerateRAM=50000000000

# zip genome and .gtf
bgzip -@ $THREADS ${IN_GENOME%.gz}
pigz -p $THREADS ${IN_GTF%.gz}

