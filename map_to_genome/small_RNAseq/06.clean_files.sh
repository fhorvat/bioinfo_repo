#/bin/bash
# ----------------Commands------------------- #
### create dirs
mkdir \
0_scripts \
1_logs \
2_original_mapping \
3_filtered_reads \
4_library_size \
5_class_reads \
6_tracks

### move files
# 0. scripts
mv \
0*.sh \
0*.R \
0_scripts/

# 1. logs
mv \
*.out \
*.tab \
Rscript*.R \
pbs*.* \
1_logs

# 2. original mapping
find . -maxdepth 1 -name "*.bam*" -not -name "*21to23nt*" -exec mv {} 2_original_mapping/ \;
mv *Unmapped.out.mate1 2_original_mapping

# 3. filtered reads
#find . -maxdepth 1 -name "*.21to23nt.bam*" -exec mv {} 3_filtered_reads/ \;

# 4. library sizes
mv \
library_hist* \
library_sizes.txt \
4_library_size

# 5. read classes
mv read_class*txt 5_class_reads

# 6. tracks
mv *.bw 6_tracks
