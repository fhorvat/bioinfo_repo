#!/bin/bash
# gets number of reads in all .fastq files in current dir
for fastq in `ls *.fastq` 
do
echo $fastq `awk '{s++}END{print s/4}' $fastq`
done
