#/bin/bash
# ----------------Loading variables------------------- #
# scripts to change 
SCRIPTS=(`ls 0*.sh 0*.R | grep -v "00b.change_scripts.sh"`)

# variables to change in scripts
OUT_PATH=`echo $PWD`

BASE_DIR=`basename ${PWD}`
GENOME_DIR=../../genomes/${BASE_DIR}
SJDB_OVERHANG=sjdbOverhang_249

INPUT_DIR=../../Raw/${BASE_DIR}/Cleaned/trimmed
HEX_PATH=~/public_html/Svoboda/bw_tracks/maternal_transcriptomes/${BASE_DIR}
LOBSANG_PATH=${PWD/common/common-lobsang}

# get number of samples
SAMPLE_LIST=(`ls $INPUT_DIR/*txt.gz | grep -v "all"`)
N_SAMPLES=`printf "%s\n" "${SAMPLE_LIST[@]%_[1-2].txt.gz}" | sort -u | wc -l | awk '{print $1-1}'`

# ----------------Commands------------------- #
# change scripts in loop
for script in "${SCRIPTS[@]}"
do 
	echo $script
	perl -pi -e "s|%OUT_PATH|$OUT_PATH|" $script
	perl -pi -e "s|%N_SAMPLES|$N_SAMPLES|" $script
	perl -pi -e "s|%GENOME_DIR|$GENOME_DIR|" $script
	perl -pi -e "s|%SJDB_OVERHANG|$SJDB_OVERHANG|" $script
	perl -pi -e "s|%BASE_DIR|$BASE_DIR|" $script
	perl -pi -e "s|%HEX_PATH|$HEX_PATH|" $script
	perl -pi -e "s|%LOBSANG_PATH|$LOBSANG_PATH|" $script
done
