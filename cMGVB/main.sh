#!/usr/bin/bash

# make
./make_cu.sh

# Move to data
cd ../data

# Run mgvb_cuda
time ../code/mgvb_2026_cuda ECL08730.raw.ms2.dso.mms sequences1.db config_focused.rms 99276 1

# Filter py PEP
Rscript ../code/PEP_for_gMGVB.R

# Copy results.txt to results
cp sig_results.txt ../results/. 

# Process in sqlite3
sqlite3 results.db ".read ../code/sqlite_script.sql"
cp modified_peptides.txt ../results/.
rm results.db