#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.fido_download
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set URL
URL="http://fido.iabio.eu/xxx"

# ----------------Commands------------------- #
# download data
wget --show-progress --recursive --no-host-directories --reject="index.html" ${URL}/
 
