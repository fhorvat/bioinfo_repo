#/bin/bash
# ----------------Loading variables------------------- #
#### experiment
EXPERIMENT=${PWD%/Data*}
EXPERIMENT=${EXPERIMENT##/*/}
EXPERIMENT_NAME="testis.8.5dpp"

#### input
SINGLE_END=TRUE
IN_GENOME=golden_hamster/MesAur1.0.GCA_000349665.1
INPUT_DIR=../../Raw/Links
ENSEMBL_VERSION=99
SJDB_OVERHANG=sjdbOverhang_100
HEX_OUT=in_house/hamster_KO/${EXPERIMENT}

# input fastq files
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.txt.gz" -not -name "*all*" \)))

### load
# scripts to change
SCRIPTS=($(find . -maxdepth 1 \( -name "0*.sh" -not -name "00b.change_scripts.sh" \)))

# input genome files
GENOME_BASE=/common/DB/genome_reference
GENOME_PATH=${GENOME_BASE}/${IN_GENOME}
STAR_INDEX=${GENOME_PATH}/STAR_index.2.7/${SJDB_OVERHANG}
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

# input features for classifying reads
FEATURES_EXONS=$(find ${GENOME_PATH} -maxdepth 1 -name "ensembl.${ENSEMBL_VERSION}.*.UCSCseqnames.reducedExons.RDS")
FEATURES_RMSK=$(find $GENOME_PATH -maxdepth 1 -name "rmsk*clean.fa.out.gz")

# paths for tracks
HEX_PATH=~/public_html/Svoboda/bw_tracks/${HEX_OUT}
LOBSANG_PATH=${PWD/common/common-lobsang}
