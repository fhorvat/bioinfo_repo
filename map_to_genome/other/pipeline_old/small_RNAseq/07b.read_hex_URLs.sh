#/bin/bash
# ----------------Loading variables------------------- #
EXPERIMENT=Yang_2016_SciAdv_GSE83581
HEX_PATH=~/public_html/Svoboda/bw_tracks/accessory_data_sets/${EXPERIMENT}
TABLE_PATH=http://hex.bioinfo.hr/~fhorvat/${HEX_PATH#~/public_html}/log.tracks_URL.txt

# ----------------Commands------------------- #
# get table from hex
wget $TABLE_PATH

# change permissions
find ./6_tracks -name "*.bw" -exec chmod 744 {} \;
find . -name "*.bam*" -exec chmod 744 {} \;
