#/bin/bash
# ----------------Loading variables------------------- #
# scripts to change 
SCRIPTS=($(find . -name "0*.sh" -or -name "0*.R" | grep -v "00b.change_scripts.sh"))

# variables to change in scripts
EXPERIMENT="CNOT6L"
EXPERIMENT_NAME="CNOT6L"

GENOME_PATH=/common/DB/genome_reference
GENOME_DIR=${GENOME_PATH}/mouse/mm10.GRCm38.GCA_000001635.2
SJDB_OVERHANG=sjdbOverhang_249

# paths for tracks
HEX_PATH=~/public_html/Svoboda/bw_tracks/in_house/RNAi_piRNA_paper/${EXPERIMENT}/perfect_alignments.all_multimappers
LOBSANG_PATH=${PWD/common/common-lobsang}

# out path
OUT_PATH=`echo $PWD`

# get number of samples
SAMPLE_LIST=($(find ../../Raw/Links -maxdepth 1 -name "*.txt.gz"))
N_SAMPLES=($(printf "%s\n" "${SAMPLE_LIST[@]%_[1-2].txt.gz}" | sort -u | wc -l | awk '{print $1-1}'))

# ----------------Commands------------------- #
# change scripts in loop
for script in "${SCRIPTS[@]}"
do 
	echo $script
	chmod 744 $script
	perl -pi -e "s|%OUT_PATH|$OUT_PATH|" $script
	perl -pi -e "s|%N_SAMPLES|$N_SAMPLES|" $script
	perl -pi -e "s|%GENOME_DIR|$GENOME_DIR|" $script
	perl -pi -e "s|%SJDB_OVERHANG|$SJDB_OVERHANG|" $script
	perl -pi -e "s|%HEX_PATH|$HEX_PATH|" $script
	perl -pi -e "s|%LOBSANG_PATH|$LOBSANG_PATH|" $script
	perl -pi -e "s|%EXPERIMENT_NAME|$EXPERIMENT_NAME|" $script
	perl -pi -e "s|%EXPERIMENT|$EXPERIMENT|" $script
done

# create soft link to the bbmap reference
find ${GENOME_DIR}/bbmap_index -name "ref" -type d -exec ln -s {} \;
