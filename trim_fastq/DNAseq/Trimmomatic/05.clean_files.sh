#!/bin/bash

# ----------------Loading variables------------------- #
IN_DIR=.
IN_SEQ=(`find ${IN_DIR} -name "*.txt.gz"`)

# rename and delete files 
for FILE in ${IN_SEQ[@]}; do
	mv $FILE ${FILE/.trim/}
done
