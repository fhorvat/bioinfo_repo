#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.align.minimap2.PE
#PBS -j oe
#PBS -J 0-3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=20G

IN_DIR=../../../../../../Raw/Cleaned
IN_SEQ=(${IN_DIR}/s_mesAur_fragment*.PE_1.txt.gz)
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*.txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

IN_FASTA=($(find . -name "*fa"))

INSERT_LOW=50
INSERT_HIGH=300

MINIMAP2_PAR="-x sr \
-a \
-F$INSERT_HIGH \
--sam-hit-only \
--secondary=no \
-t $THREADS \
-k19 -w10 -A3 -B5 -O6,26"

# ----------------Commands------------------- #
# align reads to transcriptome
minimap2 \
$MINIMAP2_PAR \
$IN_FASTA \
${FILE}_1.txt.gz \
${FILE}_2.txt.gz | \
samtools sort -@ $THREADS -o ${BASE}.bam

# Short single-end reads without splicing
# -k21 # k-mer length
# -w11 # windows size
# --sr # short-read alignment heuristics
# --frag=yes # enable the fragment mode
# -A2 # matching score
# -B8 # mismatching penalty
# -O12,32 # gap open penalty
# -E2,1 # gap extension penalty
# -r50 # Bandwidth used in chaining and DP-based alignment [500]. This option approximately controls the maximum gap size. 
# -p.5 # Minimal secondary-to-primary score ratio to output secondary mappings
# -N20 # Output at most INT secondary alignments
# -f1000,5000 # ignore minimizers occuring more than INT1 times. INT2 is only effective in the --sr or -xsr mode, which sets the threshold for a second round of seeding
# -n2 # Discard chains consisting of N number of minimizers
# -m20 # Discard chains with chaining score < N
# -s40 # Minimal peak DP alignment score to output
# -g200 # Stop chain enlongation if there are no minimizers within INT-bp [10000]
# -K50m # Number of bases loaded into memory to process in a mini-batch
# --heap-sort=yes # If yes, sort anchors with heap merge, instead of radix sort
# --secondary=no # Whether to output secondary alignments
