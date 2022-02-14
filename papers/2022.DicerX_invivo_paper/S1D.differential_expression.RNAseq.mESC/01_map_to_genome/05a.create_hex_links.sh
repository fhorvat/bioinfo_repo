#/bin/bah

# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/in_house/DicerX_paper/datasets/mESC_DicerX.2018_June.RNAseq
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/DicerX_paper/datasets/mESC_DicerX.2018_June.RNAseq/Data/Mapped/STAR_mm10

# ----------------Commands------------------- #
# create links
mkdir -p ${HEX_PATH}
cd ${HEX_PATH}
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bw" \) -exec ln -s {} $HEX_PATH \;
find ${LOBSANG_PATH} -maxdepth 1 \( -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" \) -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
