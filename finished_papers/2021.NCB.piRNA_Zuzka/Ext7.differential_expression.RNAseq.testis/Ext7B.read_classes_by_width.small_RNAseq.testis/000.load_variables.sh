#/bash
# ----------------Loading variables------------------- #
#### experiment
EXPERIMENT=${PWD%/Data*}
EXPERIMENT=${EXPERIMENT##/*/}
EXPERIMENT_NAME=${EXPERIMENT}

#### input
SINGLE_END=TRUE
IN_GENOME=golden_hamster/MesAur1.0.GCA_000349665.1
INPUT_DIR=./bam_files

# input files
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.bam" \)))

### load
# scripts to change
SCRIPTS=($(find . -maxdepth 1 \( -name "0*.sh" -not -name "00b.change_scripts.sh" \)))

# input genome files
GENOME_BASE=/common/DB/genome_reference
GENOME_PATH=${GENOME_BASE}/${IN_GENOME}

# input features for classifying reads
FEATURES_RMSK=$(find $GENOME_PATH -maxdepth 1 -name "rmsk*clean.fa.out.gz")
FEATURES_ENSEMBL=$(find ${GENOME_PATH} -maxdepth 1 -name "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf")
