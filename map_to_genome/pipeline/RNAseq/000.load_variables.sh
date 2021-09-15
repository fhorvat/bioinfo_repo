#/bin/bash
# ----------------Loading variables------------------- #
#### experiment
EXPERIMENT=${PWD%/Data*}
EXPERIMENT=${EXPERIMENT##/*/}
EXPERIMENT_NAME=${EXPERIMENT%_*_*}

#### input
SINGLE_END=TRUE
IN_GENOME=mouse/mm10.GRCm38.GCA_000001635.2
INPUT_DIR=../../Raw/Links
ENSEMBL_VERSION=99
SJDB_OVERHANG=sjdbOverhang_100
HEX_OUT=accessory_data_sets/${EXPERIMENT}

# input fastq files
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.txt.gz" \)))

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
