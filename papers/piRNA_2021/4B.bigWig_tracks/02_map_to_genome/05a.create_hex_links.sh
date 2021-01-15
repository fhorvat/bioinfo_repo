#/bin/bah

# ----------------Loading variables------------------- #
EXPERIMENT=hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq
PERFECT=TRUE

# hex path
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/in_house/hamster_KO/Siomi
HEX_PATH=${HEX_PATH}/${EXPERIMENT}/STAR_multimappers

if [ ${PERFECT} == TRUE ]
then
   HEX_PATH=${HEX_PATH}/merged_replicates/perfect_reads
else
   HEX_PATH=${HEX_PATH}/merged_replicates
fi

# lobsang path
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets
LOBSANG_PATH=${LOBSANG_PATH}/${EXPERIMENT}/Data/Mapped/STAR_Siomi.multimappers

if [ ${PERFECT} == TRUE ]
then
   LOBSANG_PATH=${LOBSANG_PATH}/5_merged_replicates/perfect_reads
else
   LOBSANG_PATH=${LOBSANG_PATH}/5_merged_replicates
fi

echo $HEX_PATH

# ----------------Commands------------------- #
# create links
mkdir -p ${HEX_PATH}
cd ${HEX_PATH}
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bw" \) -exec ln -s {} $HEX_PATH \;
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" \) -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
