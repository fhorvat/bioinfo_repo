#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.aspera.EMBL
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
ASP_EXE="/common/WORK/fhorvat/programi/Aspera/cli/bin/aspera"
DATA_NAME="Malik_CCRVDANXX"
ASP_PAR="--user=malikr@img.cas.cz \
--host=faspex.embl.de \
--password=Kralik74! \
--insecure"

# ----------------Commands------------------- #
# get Faspex URL
DATA_URL=`$ASP_EXE faspex list $ASP_PAR --inbox | grep "^Faspex URL: " | grep $DATA_NAME | awk '{print $NF}'`

# download data
$ASP_EXE faspex get $ASP_PAR --file="." --url=$DATA_URL
