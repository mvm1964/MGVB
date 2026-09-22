#!/usr/bin/bash

# make
#./make_cu.sh

# Move to data
#cd ../data

# Run mgvb_cuda
time ./mgvb_2026_cuda  $1 $2 config_focused.rms $3 $4

# Filter py PEP
Rscript PEP_for_gMGVB.R

# Copy results.txt to results
#cp sig_results.txt ../results/.

# copy $2 to temp db
cp $2 sequences_tmp.db

# Process in sqlite3
sqlite3 results.db ".read sqlite_script.sql"
#cp modified_peptides.txt ../results/.

# cleanup
rm results.db
rm sequences_tmp.db
