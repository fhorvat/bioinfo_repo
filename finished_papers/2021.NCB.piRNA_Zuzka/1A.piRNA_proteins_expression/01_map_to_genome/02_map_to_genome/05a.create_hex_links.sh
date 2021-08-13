#/bin/bah

# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/accessory_data_sets/Hendrickson_2017_NatGenet_GSE72379
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Hendrickson_2017_NatGenet_GSE72379/Data/Mapped/STAR_hg38

# ----------------Commands------------------- #
# create links
mkdir -p ${HEX_PATH}
cd ${HEX_PATH}
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bw" \) -exec ln -s {} $HEX_PATH \;
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" \) -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
