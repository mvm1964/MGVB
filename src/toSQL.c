/*
 * A library of functions to populate sequences.db 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sqlite3.h"
#include "pepMass.h"

int toSQL_proteins(sqlite3 *db, char *data_tbl) {
    char str[4096];
    char data_str[150000];

    //sqlite3 *db;
    FILE *d_fl;    
    sqlite3_stmt *ppStmt;
    
    // declare vars for insert statements
    int prot_id;
    char gene[1024];
    char header[2048];
    char sequence[100000];
    int decoy;
    int contam; 

    // create table proteins    
    sprintf(str, "create table if not exists proteins(prot_id int, gene text, header text, sequence text, decoy int, contam int);");
    sqlite3_prepare_v2(db, str, -1, &ppStmt, 0); 
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);
    
    // open data file for reading
    if ((d_fl = fopen(data_tbl, "r")) == NULL) {
        printf("Error opening data file!\n");
        exit(1);
    }

    if (sqlite3_prepare_v2(db, "insert into proteins values(?, ?, ?, ?, ?, ?)", -1, &ppStmt, NULL)) {
        printf("Error executing sql statement\n");
        sqlite3_close(db);
        exit(-1);
    }

    char tokens[10][100000];

    // will use getline()
    //char *line = NULL;

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    while (fgets(data_str, sizeof(data_str), d_fl) != NULL) {
        //sscanf(data_str, "%d\t%s\t%s\t%s\t%d\t%d", &prot_id, gene, header, sequence, &decoy, &contam);
        // sscanf does not work well with strings with spaces, should use strtok instead
    //while (getline(&line, &len, fp) != -1) {
        char *token = strtok(data_str, "\t");
        //prot_id = atoi(token);
        int i = 0;
        while (token != NULL) {
            //printf("%s\n", token);
            strcpy(tokens[i++], token);
            token = strtok(NULL, "\t");
        }
        //exit(1);
        prot_id = atoi(tokens[0]);
        strcpy(gene, tokens[1]);
        strcpy(header, tokens[2]);
        strcpy(sequence, tokens[3]);
        decoy = atoi(tokens[4]);
        contam = atoi(tokens[5]);
                
        sqlite3_bind_int(ppStmt, 1, prot_id);
        sqlite3_bind_text(ppStmt, 2, gene, -1, NULL);
        sqlite3_bind_text(ppStmt, 3, header, -1, NULL);
        sqlite3_bind_text(ppStmt, 4, sequence, -1, NULL);
        sqlite3_bind_int(ppStmt, 5, decoy);
        sqlite3_bind_int(ppStmt, 6, contam);
        sqlite3_step(ppStmt);
        sqlite3_reset(ppStmt);
    }
    sqlite3_finalize(ppStmt);    
    sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, NULL);
    
    //sprintf(str, "sqlite3 %s 'create table if not exists proteins(prot_id int, gene text, header text, sequence text, decoy int, contam int);' '.exit'", argv[1]);
    //printf("%s\n", str);
    //system(str);
    //sprintf(str, "sqlite3 %s '.separator \"\\t\"' '.import %s proteins' '.exit'", argv[1], argv[2]);
    //printf("%s\n", str);
    //system(str);

    // Update 14.09.2020: create index for post-processing
    //sprintf(str, "sqlite3 %s 'create index prot_id on proteins(prot_id, gene, header);' '.exit'", argv[1]);
    //system(str);
    if (sqlite3_prepare_v2(db, "create index prot_id on proteins(prot_id);", -1, &ppStmt, NULL)) {
        printf("Error preparing sql statement\n");
        sqlite3_close(db);
        exit(-1);
    }
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);
    fclose(d_fl);

    return 0;
}

int toSQL_peptides(sqlite3 *db, char *data_tbl) {
    char str[4096];
    char data_str[150000];

    FILE *d_fl;
    sqlite3_stmt *ppStmt;

    // declare vars for insert statements
    char sequence[100000];
    double mass;
    int prot_id;
    int missed;
    int decoy;
    int contam;

    // create table proteins
    sprintf(str, "create table if not exists raw_peptides(sequence text, mass double,  prot_id int, missed int, decoy int, contam int);");
    sqlite3_prepare_v2(db, str, -1, &ppStmt, 0);
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);
    // open data file for reading
    if ((d_fl = fopen(data_tbl, "r")) == NULL) {
        printf("Error opening data file!\n");
        exit(1);
    }

    if (sqlite3_prepare_v2(db, "insert into raw_peptides values(?, ?, ?, ?, ?, ?)", -1, &ppStmt, NULL)) {
        printf("Error executing sql statement\n");
        sqlite3_close(db);
        exit(-1);
    }

    char tokens[10][100000];

    // will use getline()
    //char *line = NULL;

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    while (fgets(data_str, sizeof(data_str), d_fl) != NULL) {
        //sscanf(data_str, "%d\t%s\t%s\t%s\t%d\t%d", &prot_id, gene, header, sequence, &decoy, &contam);
        // sscanf does not work well with strings with spaces, should use strtok instead
    //while (getline(&line, &len, fp) != -1) {
        char *token = strtok(data_str, "\t");
        //prot_id = atoi(token);
        int i = 0;
        while (token != NULL) {
            //printf("%s\n", token);
            strcpy(tokens[i++], token);
            token = strtok(NULL, "\t");
        }

        //exit(1);
        strcpy(sequence, tokens[0]);
        mass = atof(tokens[1]);
        prot_id = atoi(tokens[2]);
        missed = atoi(tokens[3]);
        decoy = atoi(tokens[4]);
        contam = atoi(tokens[5]);

        sqlite3_bind_text(ppStmt, 1, sequence, -1, NULL);
        sqlite3_bind_double(ppStmt, 2, mass);
        sqlite3_bind_int(ppStmt, 3, prot_id);
        sqlite3_bind_int(ppStmt, 4, missed);
        sqlite3_bind_int(ppStmt, 5, decoy);
        sqlite3_bind_int(ppStmt, 6, contam);
        sqlite3_step(ppStmt);
        sqlite3_reset(ppStmt);
    }
    sqlite3_finalize(ppStmt);
    sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, NULL);

    if (sqlite3_prepare_v2(db, "create table if not exists peptides as select sequence, max(mass) as 'mass', group_concat(distinct prot_id) as prot_id, max(missed) as 'missed', max(decoy) as 'decoy', max(contam) as 'contam' from raw_peptides group by sequence;", -1, &ppStmt, NULL)) {
        printf("Error preparing sql statement\n");
        sqlite3_close(db);
        exit(-1);
    }

    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);


    /*if (sqlite3_prepare_v2(db, "create index if not exists mass on peptides(mass, sequence, prot_id, missed, decoy, contam);", -1, &ppStmt, NULL)) {
        printf("Error preparing sql statement\n");
        sqlite3_close(db);
        exit(-1);
    }
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);*/

    sprintf(str, "drop table raw_peptides;");
    if (sqlite3_prepare_v2(db, str, -1, &ppStmt, NULL)) {
        printf("Error preparing sql statement to drop raw_peptides table!\n");
        sqlite3_close(db);
        exit(-1);
    }
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt); 

    fclose(d_fl);

    return 0;
}


int toSQL_mod_pep(sqlite3 *db, char *data_tbl, char *cnfg_file) {
    
    sqlite3_stmt *ppStmt;
    FILE *cnfg;
    FILE *d_fl;

    // open config file for reading
    if ((cnfg = fopen(cnfg_file, "r")) == NULL) {
        printf("Error openning config file!\n");
        exit(1);
    }

    // open data file for reading
    if ((d_fl = fopen(data_tbl, "r")) == NULL) {
        printf("Error opening data file!\n");
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

    // create table mod_pep
    sprintf(str, "create table if not exists mod_pep(pep_id int, mass double);");
    sqlite3_prepare_v2(db, str, -1, &ppStmt, 0);
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    // need to add ptm names as columns
    int j = 0;
    while (j < i) {
        sprintf(str, "alter table mod_pep add column %s int;", ptm[j].name);
        sqlite3_prepare_v2(db, str, -1, &ppStmt, 0);
        sqlite3_step(ppStmt);
        sqlite3_finalize(ppStmt);
        j++;
    }
    strcpy(str, "insert into mod_pep values(?, ?");
    j = 0;
    while (j < i) {
        strcat(str, ", ?");
        j++;
    }     
    strcat(str, ")");

    // for debugging
    //printf("%s\n", str);

    if (sqlite3_prepare_v2(db, str, -1, &ppStmt, NULL)) {
        printf("Error preparing sql statement for mod_pep loading to sqlite db\n");
        sqlite3_close(db);
        exit(-1);
    }
    int pep_id;
    double mass;
    char tokens[2 + i][1024];
    char data_str[1024];
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    while (fgets(data_str, sizeof(data_str), d_fl) != NULL) {
        char *token = strtok(data_str, "\t");
        int t = 0;
        while (token != NULL) {
            strcpy(tokens[t++], token);
            token = strtok(NULL, "\t");
        }
        pep_id = atoi(tokens[0]);
        mass = atof(tokens[1]);
        sqlite3_bind_int(ppStmt, 1, pep_id);
        sqlite3_bind_double(ppStmt, 2, mass);
        for (int w = 0; w < i; w++) {
            sqlite3_bind_int(ppStmt, 3 + w, atoi(tokens[2 + w]));
        }
        sqlite3_step(ppStmt);
        sqlite3_reset(ppStmt);
    }
    sqlite3_finalize(ppStmt);
    sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, NULL);

    sprintf(str, "create index if not exists mod_mass on mod_pep(mass);");
    sqlite3_prepare_v2(db, str, -1, &ppStmt, 0);
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    sprintf(str, "vacuum;");
    sqlite3_prepare_v2(db, str, -1, &ppStmt, 0);
    sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    fclose(d_fl);
    fclose(cnfg);
    free(ptm);
    
    return 0;
}



