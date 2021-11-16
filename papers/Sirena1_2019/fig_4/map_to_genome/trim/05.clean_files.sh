#!/bin/bash

# ----------------Loading variables------------------- #
IN_DIR=.
IN_SEQ=(`ls ${IN_DIR}/*.clean.trim.txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*.clean.trim.txt.gz}" | sort -u`)

# rename and delete files 
for FILE in ${UNIQ_SEQ[@]}; do
	BASE=${FILE#${IN_DIR}/}	
	rm ${BASE}_1.trim.txt.gz ${BASE}_2.trim.txt.gz
	rm ${BASE}_1.dirty.trim.txt.gz ${BASE}_2.dirty.trim.txt.gz
	mv ${BASE}_1.clean.trim.txt.gz ${BASE}_1.txt.gz
	mv ${BASE}_2.clean.trim.txt.gz ${BASE}_2.txt.gz
done

# add links
IN_SEQ=(`ls ${IN_DIR}/*.txt.gz`)
for FILE in ${IN_SEQ[@]};do
	ln -s $FILE ${FILE/txt/fq}
done

### create dirs
mkdir \
0_scripts \
3_logs

### move files
# original scripts
mv \
0*.sh \
0*.R \
0*.txt \
0_scripts/

# other logs and scripts
mv \
*.log \
*.stats \
*.fasta \
Rscript*.R \
pbs*.* \
3_logs/
