#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05.full_max_likehood_reconstruction.RAxML
#PBS -l select=ncpus=6:mem=1g
##PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input 
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.fasta" -or -name "*.fasta-gb" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}

# get files
FILE=Dicer1.Metazoa.OrthoDB.20211117.taxonomy_and_sequences_PS_CDD.ClustalOmega.msa.fasta-gb
FILE_TREE=${FILE}.tree
FILE_MODEL=${FILE}.AIC.best_fit.txt

# get base
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE}

# script to run
SCRIPT=/common/WORK/fhorvat/programi/RAxML/standard-RAxML/raxmlHPC-PTHREADS-SSE3

# ----------------Commands------------------- #
### get best-fit models of protein evolution from prottest3
# result: 
# "^[A-Z]+(?=\\+)" is amino acid replacement matrix (one of the JTT, LG, DCMut, MtREV, MtMam, MtArt, Dayhoff, WAG, RtREV, CpREV, Blosum62, VT, HIVb, HIVw, FLU)
# +I models with a proportion of invariable sites
# +G models with rate variation among sites and number of categories
# +F models with empirical frequency estimation
MODEL=$(grep "Best model according to AIC:" ${FILE_MODEL} | sed 's/^Best model according to AIC: //')

### get final model of protein substitution parameter for RAxML
# Available AA substitution models: 
# DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG,  
# MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, STMTREV, DUMMY, DUMMY2, AUTO, 
# LG4M, LG4X, PROT_FILE, GTR_UNLINKED, GTR 
# With the optional "F" appendix you can specify if you want to use empirical base frequencies.

# get AA matrix - get uppercase model and remove everything after first "+"
AA_MATRIX=$(echo ${MODEL^^} | sed 's/+.*//')

# paste together
AA_SUB=$(echo "PROTGAMMA""I"${AA_MATRIX})

# consensus tree
${SCRIPT} -m ${AA_SUB}
