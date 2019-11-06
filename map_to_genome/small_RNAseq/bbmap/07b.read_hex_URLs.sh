#/bin/bash
# ----------------Loading variables------------------- #
HEX_PATH=~/public_html/Svoboda/bw_tracks/in_house/hamster_KO/hamster_testis_Mov10l.small_RNAseq/perfect_alignments.all_multimappers
TABLE_PATH=http://hex.bioinfo.hr/~fhorvat/${HEX_PATH#~/public_html}/log.tracks_URL.txt

# ----------------Commands------------------- #
# get table from hex
wget $TABLE_PATH

# change permissions
find ./6_tracks -name "*.bw" -exec chmod 744 {} \;
find . -name "*.bam*" -exec chmod 744 {} \;
