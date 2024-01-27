
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sqlite3.h"
#include <libgen.h>
#include "pepMass.h"


typedef struct {
    char name[1024];
    char experiment[1024];
} DB_FILES;	

// computes pep
int compute_PEP(sqlite3 *db, char *cond) {

    char sql_inst[4096];
    
    // compute means and sd for forward sequences
    double m_s_f, m_l_f, sd_s_f, sd_l_f, cov_f, cov_r;
    sprintf(sql_inst, "select avg(Score_1), avg(length(Sequence)), \
		    sqrt(AVG((Score_1 - 50)*(Score_1 - 50)) - AVG(Score_1 - 50)*AVG(Score_1 - 50)), \
		    sqrt(AVG((length(Sequence) - 12)*(length(Sequence) - 12)) - AVG(length(Sequence) - 12)*AVG(length(Sequence) - 12)), \
		    AVG((Score_1 - 50)*(length(Sequence) - 12)) - AVG(Score_1 - 50)*AVG(length(Sequence) - 12) \
		    from merged_aggregated where Decoy = 0 and %s;", cond);

    // debug
    //printf("%s\n", sql_inst);
    sqlite3_stmt *ppStmt;
    int rc;

    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    while ((rc = sqlite3_step(ppStmt)) != 101) {
        if (rc == 100) {
            m_s_f = sqlite3_column_double(ppStmt, 0);
            m_l_f = sqlite3_column_double(ppStmt, 1);
            sd_s_f = sqlite3_column_double(ppStmt, 2);
            sd_l_f = sqlite3_column_double(ppStmt, 3);
	    cov_f = sqlite3_column_double(ppStmt, 4);
        }
    }
    sqlite3_finalize(ppStmt);

    // debug
    //printf("%lf, %lf, %lf, %lf\n", m_s_f, m_l_f, sd_s_f, sd_l_f);

    // compute means and sd for reverse sequences
    double m_s_r, m_l_r, sd_s_r, sd_l_r;
    sprintf(sql_inst, "select avg(Score_1), avg(length(Sequence)), \
		    sqrt(AVG((Score_1 - 50)*(Score_1 - 50)) - AVG(Score_1 - 50)*AVG(Score_1 - 50)), \
		    sqrt(AVG((length(Sequence) - 12)*(length(Sequence) - 12)) - AVG(length(Sequence) - 12)*AVG(length(Sequence) - 12)), \
		    AVG((Score_1 - 50)*(length(Sequence) - 12)) - AVG(Score_1 - 50)*AVG(length(Sequence) - 12)\
		    from merged_aggregated where Decoy = 1 and %s;", cond);
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    while ((rc = sqlite3_step(ppStmt)) != 101) {
        if (rc == 100) {
            m_s_r = sqlite3_column_double(ppStmt, 0);
            m_l_r = sqlite3_column_double(ppStmt, 1);
            sd_s_r = sqlite3_column_double(ppStmt, 2);
            sd_l_r = sqlite3_column_double(ppStmt, 3);
	    cov_r = sqlite3_column_double(ppStmt, 4);
        }

    }
    sqlite3_finalize(ppStmt);
    double rho_f = cov_f/(sd_s_f*sd_l_f);
    double rho_r = cov_r/(sd_s_r*sd_l_r);

    // need to test if rho is nan
    if (isnan(rho_r)) rho_r = 0;
    if (isnan(rho_f)) rho_f = 0;

    // debug
    //printf("%lf, %lf, %lf, %lf\n", m_s_r, m_l_r, sd_s_r, sd_l_r);


    //sprintf(sql_inst, "UPDATE merged_aggregated SET PEP = 0.5*(exp(-( (power(Score_1 - %lf, 2))/(2*power(%lf, 2))) + ((power(length(Sequence) - %lf, 2))/(2*power(%lf, 2))) )) / ((exp(-( (power(Score_1 - %lf, 2))/(2*power(%lf, 2))) + ((power(length(Sequence) - %lf, 2))/(2*power(%lf, 2))) )) + exp(-( (power(Score_1 - %lf, 2))/(2*power(%lf, 2))) + ((power(length(Sequence) - %lf, 2))/(2*power(%lf, 2)))) )  where %s;", m_s_r, sd_s_r, m_l_r, sd_l_r, m_s_r, sd_s_r, m_l_r, sd_l_r, m_s_f, sd_s_f, m_l_f, sd_l_f, cond);

    sprintf(sql_inst, "UPDATE merged_aggregated SET PEP = (exp(-(power(Score_1 - %lf, 2)/(2*power(%lf, 2)*(1 - %lf*%lf)) - (2*%lf*(Score_1 - %lf)*(length(Sequence) - %lf)/(2*%lf*%lf*(1 - %lf*%lf))) + power(length(Sequence) - %lf, 2)/(2*power(%lf, 2)*(1 - %lf*%lf))))/(2*%lf*%lf*(sqrt(1 - %lf*%lf)))) / \
		    ((exp(-(power(Score_1 - %lf, 2)/(2*power(%lf, 2)*(1 - %lf*%lf)) - (2*%lf*(Score_1 - %lf)*(length(Sequence) - %lf)/(2*%lf*%lf*(1 - %lf*%lf))) + power(length(Sequence) - %lf, 2)/(2*power(%lf, 2)*(1 - %lf*%lf))))/(2*%lf*%lf*(sqrt(1 - %lf*%lf)))) + (exp(-(power(Score_1 - %lf, 2)/(2*power(%lf, 2)*(1 - %lf*%lf)) - (2*%lf*(Score_1 - %lf)*(length(Sequence) - %lf)/(2*%lf*%lf*(1 - %lf*%lf))) + power(length(Sequence) - %lf, 2)/(2*power(%lf, 2)*(1 - %lf*%lf))))/(2*%lf*%lf*(sqrt(1 - %lf*%lf))))) where %s;",
		    m_s_r, sd_s_r, rho_r, rho_r, rho_r, m_s_r, m_l_r, sd_s_r, sd_l_r, rho_r, rho_r, m_l_r, 
		    sd_l_r, rho_r, rho_r, sd_s_r, sd_l_r, rho_r, rho_r,
		    m_s_r, sd_s_r, rho_r, rho_r, rho_r, m_s_r, m_l_r, sd_s_r, sd_l_r, rho_r, rho_r, m_l_r, 
		    sd_l_r, rho_r, rho_r, sd_s_r, sd_l_r, rho_r, rho_r,
		    m_s_f, sd_s_f, rho_f, rho_f, rho_f, m_s_f, m_l_f, sd_s_f, sd_l_f, rho_f, rho_f, m_l_f, 
		    sd_l_f, rho_f, rho_f, sd_s_f, sd_l_f, rho_f, rho_f, cond);
		    
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    if (rc) {
        printf("Prepare code is %d\n", rc);
        printf("Prepared instruction is %s\n", sql_inst);
        exit(1);
    }
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);
    return rc;
}

int main(int argc, char **argv) {
    
    /*
     * Will work like this:
     *
     * 1. Open a new db called merged_db and creates a table merged_aggregated
     * 2. In a loop, attaches *.raw.ms2.mms.txt.db files and inserts their aggregated table
     * into merged_aggregated
     * 3. Computes PEP for merged_aggregated
     * 4. Copies PEP from merged aggregated to each of the attached databases
    */
    

    // get the dir name (don't forget to free)
    //char *dirName = dirname(argv[1]);
    FILE *cnfg;
    // open config file for reading
    if ((cnfg = fopen(argv[1], "r")) == NULL) {
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
    int n_files = 0;
    
    DB_FILES *db_files = malloc(sizeof(*db_files));
    if (!db_files) {
        printf("Error allocating memory for db_files array!\n");
        exit(1);
    }

    // debug
    //strcpy(file_names[0], "WT1.raw");
    //printf("%s\n", file_names[n_files]);
    //exit(1);

    // Read modifications and file names into ptm array from cnfg file.
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
	else if (strstr(str, "_rawFile")) {
            DB_FILES *tmp_names = realloc(db_files, (1 + n_files)*sizeof(DB_FILES));
            if (!tmp_names) {
                printf("Error reallocating memory for file names array!\n");
                free(tmp_names);
                exit(1);
            }

            db_files = tmp_names;

            sscanf(str, "_rawFile\t%s\t%s", db_files[n_files].name, db_files[n_files].experiment);
	    strcat(db_files[n_files].name, ".ms2.mms.txt.db");
            n_files++;
        }
    }

    // create database and various tables
    char sql_inst[4096];
    char *databaseName = "merged_db.db";
    
    // update 10.09.2021: will replace dot commands as in toSQL.c
    sprintf(sql_inst, "CREATE TABLE IF NOT EXISTS merged_aggregated(Total INTEGER, Matched INTEGER, Int_covered NUMERIC, q NUMERIC");
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
     strcat(sql_inst, ", nMod INTEGER, Missed INTEGER, Decoy INTEGER, Contam INTEGER, n_seq ITEGER, score_threshold NUMERIC, file_name TEXT, row_id INTEGER, next_max_score NUMERIC, score_diff	NUMERIC, Score_1 NUMERIC);");

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
    
    // open data file for reading
    for (int i = 0; i < n_files; i++) {
        // attach database
        sprintf(sql_inst,  "attach '%s' as fileDB;", db_files[i].name);
        rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
        rc = sqlite3_step(ppStmt);
        sqlite3_finalize(ppStmt);
        
	// insert aggregated
	sprintf(sql_inst,  "insert into merged_aggregated select * from fileDB.aggregated;");
        rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
        rc = sqlite3_step(ppStmt);
        sqlite3_finalize(ppStmt);

	sprintf(sql_inst,  "detach '%s';", "fileDB");
        rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
        rc = sqlite3_step(ppStmt);
        sqlite3_finalize(ppStmt);
    }
    
    // now compute PEP on merged_aggregated
    
    sprintf(sql_inst, "ALTER TABLE merged_aggregated ADD PEP NUMERIC;");
    rc = sqlite3_prepare_v2(db, sql_inst, -1, &ppStmt, 0);
    rc = sqlite3_step(ppStmt);
    sqlite3_finalize(ppStmt);

    // Compute PEP and add it to the aggregated table in db
    int rc1;
    /*rc1 = compute_PEP(db, "nMod = 0 and Missed = 0");
    rc1 = compute_PEP(db, "nMod = 0 and Missed = 1");
    rc1 = compute_PEP(db, "nMod = 0 and Missed > 1");

    rc1 = compute_PEP(db, "nMod = 1 and Missed = 0");
    rc1 = compute_PEP(db, "nMod = 1 and Missed = 1");
    rc1 = compute_PEP(db, "nMod = 1 and Missed > 1");

    rc1 = compute_PEP(db, "nMod = 2 and Missed = 0");
    rc1 = compute_PEP(db, "nMod = 2 and Missed = 1");
    rc1 = compute_PEP(db, "nMod = 2 and Missed > 1");

    rc1 = compute_PEP(db, "nMod > 2 and Missed = 0");
    rc1 = compute_PEP(db, "nMod > 2 and Missed = 1");
    rc1 = compute_PEP(db, "nMod > 2 and Missed > 1");
    */
    rc1 = compute_PEP(db, "nMod = 0 and Missed = 0 and Charge = 2");
    rc1 = compute_PEP(db, "nMod = 0 and Missed = 1 and Charge = 2");
    rc1 = compute_PEP(db, "nMod = 0 and Missed > 1 and Charge = 2");

    rc1 = compute_PEP(db, "nMod = 1 and Missed = 0 and Charge = 2");
    rc1 = compute_PEP(db, "nMod = 1 and Missed = 1 and Charge = 2");
    rc1 = compute_PEP(db, "nMod = 1 and Missed > 1 and Charge = 2");

    rc1 = compute_PEP(db, "nMod = 2 and Missed = 0 and Charge = 2");
    rc1 = compute_PEP(db, "nMod = 2 and Missed = 1 and Charge = 2");
    rc1 = compute_PEP(db, "nMod = 2 and Missed > 1 and Charge = 2");

    rc1 = compute_PEP(db, "nMod > 2 and Missed = 0 and Charge = 2");
    rc1 = compute_PEP(db, "nMod > 2 and Missed = 1 and Charge = 2");
    rc1 = compute_PEP(db, "nMod > 2 and Missed > 1 and Charge = 2");

    // now charge > 2 
    rc1 = compute_PEP(db, "nMod = 0 and Missed = 0 and Charge != 2");
    rc1 = compute_PEP(db, "nMod = 0 and Missed = 1 and Charge != 2");
    rc1 = compute_PEP(db, "nMod = 0 and Missed > 1 and Charge != 2");

    rc1 = compute_PEP(db, "nMod = 1 and Missed = 0 and Charge != 2");
    rc1 = compute_PEP(db, "nMod = 1 and Missed = 1 and Charge != 2");
    rc1 = compute_PEP(db, "nMod = 1 and Missed > 1 and Charge != 2");

    rc1 = compute_PEP(db, "nMod = 2 and Missed = 0 and Charge != 2");
    rc1 = compute_PEP(db, "nMod = 2 and Missed = 1 and Charge != 2");
    rc1 = compute_PEP(db, "nMod = 2 and Missed > 1 and Charge != 2");

    rc1 = compute_PEP(db, "nMod > 2 and Missed = 0 and Charge != 2");
    rc1 = compute_PEP(db, "nMod > 2 and Missed = 1 and Charge != 2");
    rc1 = compute_PEP(db, "nMod > 2 and Missed > 1 and Charge != 2");
    
    sqlite3_close(db);     

    free(ptm);

    //for (int i = 0; i < 300; i++) free(file_names[i]);
    free(db_files);
    fclose(cnfg);
    return 0;

}   

