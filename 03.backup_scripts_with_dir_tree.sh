cd /common/WORK/fhorvat/Projekti/Svoboda
find . \( -name "*sh" -or -name "0*txt" \) -print0 | xargs -0 tar -zcf ./scripts/scripts_in_dirs.tar.gz
