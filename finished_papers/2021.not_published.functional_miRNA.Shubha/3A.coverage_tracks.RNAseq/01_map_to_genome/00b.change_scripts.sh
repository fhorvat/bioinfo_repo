#/bin/bash
# ----------------Loading variables------------------- #
# scripts to change 
SCRIPTS=(`ls 0*.sh 0*.R | grep -v "00b.change_scripts.sh"`)

# variables to change in scripts
OUT_PATH=`echo $PWD`
GENOME_PATH=/common/DB/genome_reference
GENOME_DIR=${GENOME_PATH}/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2
SJDB_OVERHANG=sjdbOverhang_249
HEX_PATH=~/public_html/Svoboda/bw_tracks/accessory_data_sets/Graf_2014_ProcNatlAcadSciUSA_GSE52415.bosTau9
LOBSANG_PATH=${PWD/common/common-lobsang}

# get number of samples
SAMPLE_LIST=(`ls ../../Raw/Links/*.txt.gz`)
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
	perl -pi -e "s|%HEX_PATH|$HEX_PATH|" $script
	perl -pi -e "s|%LOBSANG_PATH|$LOBSANG_PATH|" $script
done
