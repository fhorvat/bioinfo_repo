# dirs and paths
SAMPLE_IN=/common/DB/SRA/sra/Svoboda/2017_download/Zuzka
SAMPLE_OUT=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/Zuzka

# make list of all SRA studies
ls -l $SAMPLE_IN | grep "SRP" | awk '{print $9}' > sra_studies.txt

# read dirs from list, make dir structure, make links to raw .fastq
for SRA_IN in `cat ${SAMPLE_OUT}/sra_studies.txt`;do
echo $SRA_IN
mkdir ${SAMPLE_OUT}/${SRA_IN} ${SAMPLE_OUT}/${SRA_IN}/Data ${SAMPLE_OUT}/${SRA_IN}/Data/Raw ${SAMPLE_OUT}/${SRA_IN}/Data/Raw/QC ${SAMPLE_OUT}/${SRA_IN}/Data/Raw/Cleaned ${SAMPLE_OUT}/${SRA_IN}/Data/Raw/Links ${SAMPLE_OUT}/${SRA_IN}/Data/Mapped ${SAMPLE_OUT}/${SRA_IN}/Data/Mapped/STAR_mm10
ln -s ${SAMPLE_IN}/${SRA_IN}/*.fastq.gz ${SAMPLE_OUT}/${SRA_IN}/Data/Raw/Links
done

