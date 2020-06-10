#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.bowtie2_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=6

# set input variables
IN_GENOME=`ls ../*.fa.gz`
GENOME_NAME=`basename $IN_GENOME`
GENOME_NAME=${GENOME_NAME%.fa.gz}
GENOME_NAME=mm10

# ----------------Commands------------------- #
# unzip genome and .gtf
unpigz -p $THREADS ${IN_GENOME}

# generate index
bowtie2-build --threads $THREADS ${IN_GENOME%.gz} $GENOME_NAME 

# zip genome and .gtf
pigz -p $THREADS ${IN_GENOME%.gz}
