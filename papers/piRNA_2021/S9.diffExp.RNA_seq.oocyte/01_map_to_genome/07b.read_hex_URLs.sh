#/bin/bash
# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/in_house/hamster_KO/hamster_oocyte_Mov10l.RNAseq
TABLE_PATH=http://hex.bioinfo.hr/~fhorvat/${HEX_PATH#~/public_html}/log.tracks_URL.txt

# ----------------Commands------------------- #
# get table from hex
wget $TABLE_PATH

# change permissions
find . -name "*.bw" -exec chmod 744 {} \;
find . -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" -exec chmod 744 {} \;
