#include <stdio.h>
#include <libgen.h>
#include <stdlib.h>
#include <time.h>
//#include <omp.h>
#include <string.h>
#include "sqlite3.h"
#include "toSQL.h"

int main(int argc, char **argv) {

    sqlite3 *db;
    //sqlite3_stmt *ppStmt;
    int rc;
    rc = sqlite3_open(argv[1], &db);
    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        exit(1);
    }
    rc = toSQL_proteins(db, argv[2]);
    if(rc) {
        sqlite3_close(db);
        exit(1);
    }

    rc = toSQL_peptides(db, argv[3]);
    if(rc) {
        sqlite3_close(db);
        exit(1);
    }
    //sqlite3_close(db);
 
    char str[1024];
    sprintf(str, "mod_pep %s %s", argv[1], argv[4]);
    system(str);

    rc = toSQL_mod_pep(db, "mod_pep.txt", "config.rms");
    if(rc) {
        sqlite3_close(db);
        exit(1);
    } 
    
    // 13.12.2023: was not closing db and valgrind gave errors
    sqlite3_close(db);
       
    return 0;
}
