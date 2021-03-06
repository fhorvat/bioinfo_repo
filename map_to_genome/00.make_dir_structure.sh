#/bin/bash
# ----------------Commands------------------- #
# create dir tree
IN_DIR=.
mkdir -p \
$IN_DIR/Data/Raw/QC \
$IN_DIR/Data/Raw/Cleaned \
$IN_DIR/Data/Raw/Links \
$IN_DIR/Data/Mapped/STAR_mm10 \ 
$IN_DIR/Data/Documentation \
$IN_DIR/Analysis/expression

# copy scripts
cp /common/WORK/fhorvat/Projekti/Svoboda/scripts/map_to_genome/pipeline/RNAseq/* $IN_DIR/Data/Mapped/STAR_mm10
cp /common/WORK/fhorvat/Projekti/Svoboda/scripts/differential_expression/RNAseq/0* $IN_DIR/Analysis/expression
cp /common/WORK/fhorvat/Projekti/Svoboda/scripts/map_to_genome/QC/0* $IN_DIR/Data/Raw/QC

# set permissions
find $IN_DIR -name "0*.sh" -exec chmod 744 {} \;
find $IN_DIR -name "0*.R" -exec chmod 744 {} \;
find $IN_DIR -name "0*.txt" -exec chmod 744 {} \;
