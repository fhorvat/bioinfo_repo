#/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.bowtie2_index
#PBS -l select=ncpus=1:mem=30g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
BOWTIE2_INDEX_NAME=myco
REF=../genomes/mycoplasma_genomes.fa

# ----------------Commands------------------- #
# create index
bowtie2 bowtie2-build $REF $BOWTIE2_INDEX_NAME
