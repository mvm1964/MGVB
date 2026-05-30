create table results(specNo Number, precMass Number, Score Number, prot_rowid Number, mod_rowid Number, nmod Number, missed Number, decoy Number, pep_len Number);
.header on
.separator \t
.import sig_results.txt results
delete from results where rowid = 1;
attach database "sequences.db" as sequences;
create table tmp1 as select rowid, sequence, prot_id from sequences.peptides;
create table tmp2 as select * from tmp1 inner join results on tmp1.rowid = results.prot_rowid;
create table tmp3 as select rowid,* from sequences.mod_comb;
create table modified_peptides as select * from tmp2 inner join tmp3 on tmp2.mod_rowid = tmp3.rowid; 
.output modified_peptides.txt
select * from modified_peptides;
