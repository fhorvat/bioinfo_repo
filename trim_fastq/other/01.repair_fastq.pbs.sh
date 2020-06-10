# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.repair_fastq
#PBS -j oe
##PBS -J 0-10
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
MEMORY=50g

IN_DIR=.
IN_SEQ=($(find $IN_DIR \( -name "*.fastq" -not -name "*repaired*" \)))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*.fastq}" | sort -u))
FILE=${UNIQ_SEQ[0]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fastq}

# ----------------Commands------------------- #
repair.sh -Xmx$MEMORY in=${FILE}_1.fastq in2=${FILE}_2.fastq out=${BASE}.repaired_1.fastq out2=${BASE}.repaired_2.fastq
