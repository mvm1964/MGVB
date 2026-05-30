# Automates focussed open search
# Don't forget to copy db1_min.db
# to your working directory
# If conventional run has been done
# comment out the call to mgvb_auto.sh
# to only do focused search on the existing
# *.txt.db files.

#!/usr/bin/bash
mgvb_auto.sh
mkdir focused 
cp config_focused.rms ./focused/config.rms
cp *.txt.db ./focused/.
cp *.ms2.dso.mms ./focused/.
sqlite3 db1_min.db "attach database \"sequences.db\" as db2; create table proteins as select * from db2.proteins; create index prot_id on proteins(prot_id);"

cp db1_min.db ./focused/.
cp merged_db.db ./focused/.
cd focused

file_arr1=(*.txt.db)
file_arr2=(*.ms2.dso.mms)
for ((i=0; i < ${#file_arr1[@]} && i < ${#file_arr2[@]}; i++))
do
   deep_seq "${file_arr1[i]}" db1_min.db sequences.db proteins_from_fasta.txt peptides.txt
   nSpec=$(get_nspectra "${file_arr2[i]}")
   
   time mgvb_2026_cuda "${file_arr2[i]}" sequences.db config_focused.rms $nSpec 1

   # Filter py PEP
   Rscript PEP_for_gMGVB.R
   sqlite3 results.db ".read sqlite_script.sql"
   cp sig_results.txt "${file_arr2[i]}.txt"
   cp modified_peptides.txt "${file_arr2[i]}.mod_pep.txt" 
   rm "${file_arr1[i]}"
   rm sequences.db
done
