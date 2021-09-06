#/bin/bash
# ----------------Commands------------------- #
### create dirs
mkdir \
0_scripts \
1_mapped \
2_unmapped \
3_logs \
4_library_size

### move files
# 0. scripts
mv \
0*.sh \
0*.R \
0*.txt \
0_scripts/

# 1. mapped
#find . -maxdepth 1 \( -name "*.bam*" -not -name "*21to23nt*" -not -name "*19to32nt*" -not -name "*24to31nt*" \) -exec mv {} 2_mapped/ \;

# 2. unmapped
mv *Unmapped.out.mate1 2_unmapped

# 3. logs
mv \
*.out \
*.tab \
pbs*.* \
log.tracks_URL.txt \
log.read_stats.txt \
s_*.read_stats.txt \
3_logs

# 4. library sizes
mv \
library_hist* \
library_sizes.txt \
4_library_size
