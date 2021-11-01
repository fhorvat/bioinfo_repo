#/bin/bash
# ----------------Loading variables------------------- #
#### experiment
EXPERIMENT=${PWD%/Data*}
EXPERIMENT=${EXPERIMENT##/*/}
EXPERIMENT_NAME=${EXPERIMENT%_*_*}

#### input
SINGLE_END=TRUE
IN_GENOME=golden_hamster/Siomi_assembly.fixed
INPUT_DIR=..
SJDB_OVERHANG=sjdbOverhang_100

# input fastq files
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.txt.gz" -not -name "*atrim.txt.gz" \)))

### load
# input genome files
GENOME_BASE=/common/DB/genome_reference
GENOME_PATH=${GENOME_BASE}/${IN_GENOME}
STAR_INDEX=${GENOME_PATH}/STAR_index.2.7/${SJDB_OVERHANG}
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt
