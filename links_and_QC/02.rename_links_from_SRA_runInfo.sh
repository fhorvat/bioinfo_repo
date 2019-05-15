#!/bin/bash
# ----------------Loading variables------------------- #
TABLE=../../Documentation/*runInfo.txt
INPUT_DIR=.
IN_FASTQ=(`ls ${INPUT_DIR}/*fastq.gz`)

# ----------------Commands------------------- #
# get number of mapped reads in millions
for FILE in "${IN_FASTQ[@]}"
do

	# get base name 
	BASE=${FILE#${INPUT_DIR}/}
	BASE=${BASE%%\.*}
	BASE=${BASE%_[1-2]}

	# get extension
	EXT=${FILE#${INPUT_DIR}/}
	EXT=.${EXT#*.}
	
	# get number of file in paired-end sequencing 
	PAIR_N=${FILE#${INPUT_DIR}/}
	PAIR_N=${PAIR_N%${EXT}}
	PAIR_N=${PAIR_N#${BASE}}

	# change fastq in extension to txt
	EXT=${EXT/fastq.gz/txt.gz}
	
	# set whether sequences are paired-end or single-end
	PAIR=".PE"
	
	# go through table, grep file, set new name, remove all whitespace from new name
	NEW_NAME=`grep $BASE $TABLE | awk -F "\t" -v PAIR_N=$PAIR_N -v EXT=$EXT -v PAIR=$PAIR '{print "s_" $8 "_" $9 PAIR PAIR_N EXT}'`
	NEW_NAME=${NEW_NAME// /}
	
	# rename command
	echo "mv" $FILE $NEW_NAME >> 03.rename_fastq.sh
	
done

# set permissions
chmod 744 03.rename_fastq.sh
