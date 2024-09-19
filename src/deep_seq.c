/*
   Prepares target databases for focused scorer_mpi search

   Author: Metodi V. Metodiev
   01.03.2023
   Ravna, Bulgaria
   All rights reserved!

   Algorithm: 
   1. A focused fasta file containing only detected proteins is exported from
      the db file produced by select_by_prob: need to write a program that
      gets the proteins from db: select header, sequence from proteins where ...
      and then prints them to a file in fasta format
   2. Digest and create new sequence db in th enew directory.
    
   Need a new config.rms file with the new fasta file in it and no modifications       

*/

#include <stdio.h>
#include <libgen.h>
#include <stdlib.h>
#include <time.h>
//#include <omp.h>
#include <string.h>
#include <sqlite3.h>
#include "toSQL.h"

// Exports fasta file from sig_proteins. Calling function need to open db file
int export_fasta(sqlite3 *db, char *db2_name) {
    char sql_inst[4096];
    
    int rc;
    sqlite3_stmt *ppStmt;
    FILE *fasta_fl;

    if ((fasta_fl = fopen("focused.fasta", "w")) == NULL) {
            printf("Error opening fasta file for writing.\n");
            exit(1);
        }

    sprintf(sql_inst, "attach database \"%s\" as db2;", db2_name);
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    rc = sqlite3_finalize(ppStmt);


    sprintf(sql_inst, "select header, sequence from db2.proteins where prot_id in (select Proteins from sig_proteins);");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    while ((rc = sqlite3_step(ppStmt)) != 101) {
        if (rc == 100) {
            fprintf(fasta_fl, "%s\n", (char *)sqlite3_column_text(ppStmt, 0));
            fprintf(fasta_fl, "%s\n", (char *)sqlite3_column_text(ppStmt, 1));            
        }
        else {
            printf("Error executing sql statement: %d.\n", rc);
            exit(1);
        }


    }
    rc = sqlite3_finalize(ppStmt);
    //sqlite3_close(db);

    fclose(fasta_fl);
    return 0;
}

int main(int argc, char **argv) {

    // argv[1] is the db file generated from select_by_prob

    char sql_inst[4096];
    sqlite3_stmt *ppStmt;

    sqlite3 *db;
    //sqlite3_stmt *ppStmt;
    int rc;
    rc = sqlite3_open(argv[1], &db);
    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        exit(1);
    }

    
    rc = export_fasta(db, argv[2]);
    sqlite3_close(db);
    //exit(0);

    system("digest_universal config.rms peptides.txt");


    // now open a new db and create sequences
    rc = sqlite3_open(argv[3], &db);
    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        exit(1);
    }

    rc = toSQL_proteins(db, argv[4]);
    if(rc) {
        sqlite3_close(db);
        exit(1);
    }

    rc = toSQL_peptides(db, argv[5]);
    if(rc) {
        sqlite3_close(db);
        exit(1);
    }
    //sqlite3_close(db);
    
    /*char str[1024];
    sprintf(str, "mod_pep %s %s", argv[1], argv[4]);
    system(str);

    rc = toSQL_mod_pep(db, "mod_pep.txt", "config.rms");
    if(rc) {
        sqlite3_close(db);
        exit(1);
    }*/

    // copy mod_comb from db1_min.db
    sprintf(sql_inst, "attach database \"%s\" as db2;", argv[2]);
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    rc = sqlite3_finalize(ppStmt); 

    sprintf(sql_inst, "create table mod_comb as select * from db2.mod_comb;");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    rc = sqlite3_finalize(ppStmt);

    sprintf(sql_inst, "create index delta_mass on mod_comb(Delta_mass);");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    rc = sqlite3_finalize(ppStmt);

    sqlite3_close(db);
    return 0;
}
