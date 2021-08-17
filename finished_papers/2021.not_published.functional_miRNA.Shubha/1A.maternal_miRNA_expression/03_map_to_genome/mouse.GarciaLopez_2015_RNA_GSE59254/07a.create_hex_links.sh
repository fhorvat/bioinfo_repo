#/bin/bah

# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/accessory_data_sets/GarciaLopez_2015_RNA_GSE59254
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/GarciaLopez_2015_RNA_GSE59254/Data/Mapped/STAR_mm10

# ----------------Commands------------------- #
# create links
mkdir -p ${HEX_PATH}
cd ${HEX_PATH}
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bw" \) -exec ln -s {} $HEX_PATH \;
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" \) -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
