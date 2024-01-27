
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sqlite3.h"
#include <libgen.h>
#include "pepMass.h"


int main(int argc, char **argv) {
    
    /*  
     *  The new algorithm will work like this:
     *  1. This program will create separate db files as before but only as far as aggregated table
     *  just before PEP is computed. Will reuse compute_by_pep for this purpose.
     *  2. A new program will create a merged_aggregated db with a single table and return.
     *  3. The last program will work with the db files created in 1 aand will use the one from 2
     *  to finish the job. Will reuse the code from the rest of select_by_pep.
     *
    */ 
    
    // get the dir name (don't forget to free)
    FILE *cnfg;
    // open config file for reading
    if ((cnfg = fopen(argv[2], "r")) == NULL) {
        printf("Error openning config file!\n");
        exit(1);
    }

    // Create ptm array and populate it from config file.
    PTM *ptm = malloc(sizeof(*ptm));
    if (!ptm) {
        printf("Error allocating memory for ptm array!\n");
        exit(1);
    }

    //char *token;
    char str[2048];
    int i = 0;

    // Read modifications into ptm array from cnfg file.
    while (fgets(str, sizeof(str), cnfg) != NULL) {

        // Change code to read from new type config file
        if (strstr(str, "_varMod")) {
            PTM *tmp = realloc(ptm, (1 + i)*sizeof(PTM));
            if (!tmp) {
                printf("Error reallocating memory for ptm array!\n");
                free(ptm);
                exit(1);
            }
            ptm = tmp;

            sscanf(str, "_varMod\t%s\t%lf\t%s", ptm[i].name, &ptm[i].delta_mass, ptm[i].site);
            i++;
        }
    }

    // create database and various tables
    char sql_inst[4096];
    char databaseName[2048] = {}; // very important to initialize like this. Otherwise garbage in the name breaks the program
    strcat(databaseName, argv[1]);
    strcat(databaseName, ".db");
    
    // update 10.09.2021: will replace dot commands as in toSQL.c
    sprintf(sql_inst, "CREATE TABLE IF NOT EXISTS scans(Total INTEGER, Matched INTEGER, Int_covered NUMERIC, q NUMERIC");
    char str3[2048];
    int j = 0;
    while (j < i) {
        sprintf(str3, ", %s_probs TEXT", ptm[j].name);
        strcat(sql_inst, str3);
        j++;
    }    
    strcat(sql_inst, ", Fragments TEXT, Scan INTEGER, PrecScan INTEGER, RetTime NUMERIC, Sequence TEXT, Charge INTEGER, PrecMass NUMERIC, Score NUMERIC, Proteins TEXT");
    j = 0;
    while (j < i) {
        sprintf(str3, ", %s INTEGER", ptm[j].name);
        strcat(sql_inst, str3);
        j++;
    } 
     strcat(sql_inst, ", nMod INTEGER, Missed INTEGER, Decoy INTEGER, Contam INTEGER, n_seq ITEGER, score_threshold NUMERIC, file_name TEXT);");

    // for debugging
    printf("%s\n", sql_inst);

    // Open sqlite database
    sqlite3 *db;
    sqlite3_stmt *ppStmt;
    int rc;
    rc = sqlite3_open(databaseName, &db);
    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return(1);
    }

    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);
    
    // OK, created merged_aggregated
    
    FILE *d_fl;
    // open data file for reading
    if ((d_fl = fopen(argv[1], "r")) == NULL) {
        printf("Error opening data file!\n");
        exit(1);
    }
    
    strcpy(sql_inst, "insert into scans values(?, ?, ?, ?");
    j = 0;
    while (j < i) {
        strcat(sql_inst, ", ?");
        j++;
    }
    strcat(sql_inst,", ?, ?, ?, ?, ?, ?, ?, ?, ?");
    j = 0;
    while (j < i) {
        strcat(sql_inst, ", ?");
        j++;
    }
    strcat(sql_inst, ", ?, ?, ?, ?, ?, ?, ?)");
    // for debugging
    printf("%s\n", sql_inst);
    
    if (sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, NULL)) {
        printf("Error executing sql statement\n");
        sqlite3_close(db);
        exit(-1);
    } 

    char tokens[20 + 2*i][2048];
    char data_str[4096];
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    while (fgets(data_str, sizeof(data_str), d_fl) != NULL) {
        char *token = strtok(data_str, "\t");
        int t = 0;
        while (token != NULL) {
            strcpy(tokens[t++], token);
            token = strtok(NULL, "\t");
        }
        sqlite3_bind_int(ppStmt, 1, atoi(tokens[0]));
        sqlite3_bind_int(ppStmt, 2, atoi(tokens[1]));
	sqlite3_bind_double(ppStmt, 3, atof(tokens[2]));
	sqlite3_bind_double(ppStmt, 4, atof(tokens[3]));
        for (int w = 0; w < i; w++) {
            sqlite3_bind_text(ppStmt, 5 + w, tokens[4 + w], -1, NULL);
        }
        sqlite3_bind_text(ppStmt, 5 + i, tokens[4 + i], -1, NULL);
        sqlite3_bind_int(ppStmt, 6 + i, atoi(tokens[5 + i]));
        sqlite3_bind_int(ppStmt, 7 + i, atoi(tokens[6 + i]));
        sqlite3_bind_double(ppStmt, 8 + i, atof(tokens[7 + i]));
        sqlite3_bind_text(ppStmt, 9 + i, tokens[8 + i], -1, NULL);
        sqlite3_bind_int(ppStmt, 10 + i, atoi(tokens[9 + i])); 
        sqlite3_bind_double(ppStmt, 11 + i, atof(tokens[10 + i]));
        sqlite3_bind_double(ppStmt, 12 + i, atof(tokens[12 + i]));
        sqlite3_bind_text(ppStmt, 13 + i, tokens[13 + i], -1, NULL);
        for (int w = 0; w < i; w++) {
            sqlite3_bind_int(ppStmt, 14 + i + w, atoi(tokens[14 + i + w]));
        }
        sqlite3_bind_int(ppStmt, 14 + i + i, atoi(tokens[14 + i + i])); 
        sqlite3_bind_int(ppStmt, 15 + i + i, atoi(tokens[15 + i + i]));
        sqlite3_bind_int(ppStmt, 16 + i + i, atoi(tokens[16 + i + i]));       
        sqlite3_bind_int(ppStmt, 17 + i + i, atoi(tokens[17 + i + i]));
        sqlite3_bind_int(ppStmt, 18 + i + i, atoi(tokens[18 + i + i]));
        sqlite3_bind_double(ppStmt, 19 + i + i, atof(tokens[19 + i + i]));
        sqlite3_bind_text(ppStmt, 20 + i + i, argv[1], -1, NULL);

        sqlite3_step(ppStmt);
        sqlite3_reset(ppStmt);
    }
    sqlite3_finalize(ppStmt);
    sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, NULL);

    sqlite3_exec(db, "delete from scans where rowid = 1", NULL, NULL, NULL); 

   // need to read the config file and populate an array of ptm names
    printf("Imported %s into scans.\n", argv[1]);

    // update 17.09.2022: will delete rows with score 0
    sqlite3_exec(db, "delete from scans where Score = 0", NULL, NULL, NULL);
    printf("Deleted raws with Score = 0 from scans.\n");

    sprintf(sql_inst, "CREATE TABLE IF NOT EXISTS aggregated AS SELECT *, rowid AS row_id from scans GROUP BY Scan HAVING Score = max(Score);");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    // for debugging
    printf("Created table aggregated.\n");

    sprintf(sql_inst, "create table aggregated_2nd AS SELECT *, rowid AS row_id from scans WHERE scans.rowid NOT IN (select row_id from aggregated) GROUP BY Scan HAVING Score = max(Score);"); 
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    // for debugging
    printf("Created table aggregated_2nd.\n");

    sprintf(sql_inst, "ALTER TABLE aggregated ADD next_max_score NUMERIC;");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    // for debugging
    printf("Altered aggregated adding next_max_score.\n");
    
    sprintf(sql_inst, "create index scan on aggregated(Scan, Score);");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sprintf(sql_inst, "create index scan_2nd on aggregated_2nd(Scan, Score);");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sprintf(sql_inst, "UPDATE aggregated SET next_max_score = (select Score from aggregated_2nd where Scan = aggregated.Scan);");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sprintf(sql_inst, "UPDATE aggregated SET next_max_score = CASE WHEN next_max_score IS NOT NULL THEN next_max_score ELSE 0 END;");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    // for debugging
    printf("Updated aggregated setting next_max_score.\n");

    sprintf(sql_inst, "ALTER TABLE aggregated ADD score_diff NUMERIC;");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sprintf(sql_inst, "UPDATE aggregated SET score_diff = Score - next_max_score;"); 
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sprintf(sql_inst, "ALTER TABLE aggregated ADD Score_1 NUMERIC;");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sprintf(sql_inst, "UPDATE aggregated SET Score_1 = Score + 0.012*Score;");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sqlite3_exec(db, "delete from aggregated WHERE Score < 40 and nMod > 0", NULL, NULL, NULL);
    sqlite3_exec(db, "delete from aggregated WHERE score_diff < 6 and nMod > 0", NULL, NULL, NULL);

    sqlite3_close(db);     

    free(ptm);
    fclose(cnfg);
    fclose(d_fl);
    return 0;

}   

