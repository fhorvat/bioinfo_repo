# remove existing backup
find . -name "scripts_in_dirs.tar.gz" -exec rm {} \;

# change dir and find files, add to tar
cd /common/WORK/fhorvat/Projekti/Svoboda
find . \( -name "*sh" -or -name "0*txt" \) -exec tar -rvf ./scripts/scripts_in_dirs.tar {} \;

# gzip tar
gzip ./scripts/scripts_in_dirs.tar
