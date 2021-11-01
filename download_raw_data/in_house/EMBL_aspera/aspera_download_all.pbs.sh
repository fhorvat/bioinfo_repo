#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.aspera.EMBL
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
ASP_EXE="/common/WORK/fhorvat/programi/Aspera/cli/bin/aspera"
ASP_PAR="--user=malikr@img.cas.cz \
--host=w3-09.genecore.embl.de \
--password=Kralik74! \
--insecure"

# ----------------Commands------------------- #
# get Faspex URL
DATA_URL=(`$ASP_EXE faspex list $ASP_PAR --inbox | grep "^Faspex URL: " | awk '{print $NF}'`)

# array
DATA_URL_ONE=${DATA_URL[$PBS_ARRAY_INDEX]}

# download data
$ASP_EXE faspex get $ASP_PAR --file="." --url=$DATA_URL_ONE
