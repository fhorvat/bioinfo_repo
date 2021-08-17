#/bin/bash

# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/accessory_data_sets/Graf_2014_ProcNatlAcadSciUSA_GSE52415.bosTau9
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Graf_2014_ProcNatlAcadSciUSA_GSE52415/Data/Mapped/STAR_bosTau9

# ----------------Commands------------------- #
# create links
mkdir -p $HEX_PATH
cd $HEX_PATH
find $LOBSANG_PATH -name "*.bw" -exec ln -s {} $HEX_PATH \;
find $LOBSANG_PATH -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
