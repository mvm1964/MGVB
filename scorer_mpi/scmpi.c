#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sqlite3.h"
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "pepMass.h"
#include "item_open.h"
#include "bst2023.h"
#include "compute_site_prob_open_optimized.h"

// Struct to hold the results. Will allocate two for each spectrum to hold top 2 scorers
typedef struct Results {
    int Total;
    int Matched;
    char mod1_probs[2048];
    char mod2_probs[2048];
    char mod3_probs[2048];
    char Fragments[4096];
    char Masses[4096];
    int Scan;
    int PrecScan;
    double RetTime;
    char Sequence[1024];
    int Charge;
    double PrecMass;
    double Score;
    char Proteins[4096];
    char mod1[1024];
    char mod2[1024];
    char mod3[1024];
    int nMod;
    int Missed;
    int Decoy;
    int Contam;
    int seq_count;
    double score_threshold;
    double matched_int;
} results;

// Parser function for parsing run parameters
int parse_params(char *config_file, int *prTol, double *frTol, char fasta[5][2048], int *maxMod, int *n_losses) {

    FILE *cnfg;
    char str[1028];
    int fasta_ctr = 0;
    //int ptm_ctr = 0;

    // Open config file for reading
    if ((cnfg = fopen(config_file, "r")) == NULL) {
        printf("Error opening config file!");
        exit(1);
    }

    // Read the config file line by line and parse it.
    while (fgets(str, sizeof(str), cnfg) != NULL) {
        if (strstr(str, "_precTol")) {
            sscanf(str, "_precTol\t%d", prTol);
        }
        else if (strstr(str, "_tol")) {
            sscanf(str, "_tol\t%lf", frTol);
        }
            /*else if (strstr(str, "_q")) {
                sscanf(str, "_q\t%d", fr_q);
            }
            else if (strstr(str, "_maxSize")) {
                sscanf(str, "_maxSize\t%d", fr_maxSize);
            }*/
        else if (strstr(str, "_fasta")) {
            sscanf(str, "_fasta\t%s", fasta[fasta_ctr]); // fasta should be crated in main with sufficient size
            fasta_ctr++;
        }
        else if (strstr(str, "_maxMod")) {
            sscanf(str, "_maxMod\t%d", maxMod);
        }

	else if (strstr(str, "_nLosses")) {
            sscanf(str, "_nLosses\t%d", n_losses);
        }
        /* Having second thought about ptm. No need to rewrite so much code.
           will just modify the code in main and other functions that read ptm
        */
    }
    fclose(cnfg);
    return fasta_ctr;
}

// A file pointer for the results. Declared here so it can be used by many functions.
FILE *res;
int precTol = 10;
double tol = 0.6;
int seq_count = 0;

double get_Pscore (int n2, int n4, double q1) {   // very importnt: q1 should be double otherwise does not work
    // Variables for computing binomial coefficient and probabilities
    double res123 = 0;

    mpf_t p;
    mpf_t p1;
    mpz_t comb;
    mpf_t prob1;
    mpf_t prob2;
    mpf_t prob;
    mpf_t comb1;
    mpfr_t score;
    mpz_init_set_ui(comb, 0);
    mpf_init_set_d(prob1, 0);
    mpf_init_set_d(prob2, 0);
    mpf_init_set_d(prob, 0);

    mpf_init_set_d(p, q1/(100));
    mpf_init_set_d(p1, 1-q1/(100));
    mpf_init_set_d(comb1, 0);
    mpfr_init(score);
    mpfr_set_d(score, 0, MPFR_RNDZ);
    mpfr_t cdf;
    mpfr_init_set_d(cdf, 0, MPFR_RNDZ);

    if (n2 > 0 && n2 < n4){
        for (int i = n2; i <= n4; i++) {
            mpf_pow_ui(prob1, p, i);
            mpf_pow_ui(prob2, p1, n4-i);
            mpf_mul(prob, prob1, prob2);
            mpz_bin_uiui(comb, n4, i);
            mpf_set_z(comb1, comb);
            mpf_mul(prob, prob, comb1);
            mpfr_set_f(score, prob, MPFR_RNDZ);
            mpfr_add(cdf, cdf, score, MPFR_RNDZ);
        }

        mpfr_log10(cdf, cdf, MPFR_RNDZ);
        res123 = -10*mpfr_get_d(cdf, MPFR_RNDZ);
        // for debuging
        //printf("n2, n4, q1, sore: %d, %d, %lf, %lf\n", n2, n4, q1, res123);
    }

    else if (n2 > 0 && n2 >= n4)  {
        mpf_pow_ui(prob1, p, n4);
        mpfr_set_f(score, prob1, MPFR_RNDZ);
        mpfr_log10(score, score, MPFR_RNDZ);
        res123 = -10*mpfr_get_d(score, MPFR_RNDZ);
    }
    else res123 = 0;

    // free variables
    mpz_clear(comb);
    mpf_clear(p);
    mpf_clear(p1);
    mpf_clear(prob1);
    mpf_clear(prob2);
    mpf_clear(prob);
    mpf_clear(comb1);
    mpfr_clear(score);
    mpfr_clear(cdf);

    return res123;
}

void make_score_table(int n2, int n4, int q2, double table[][n2][n4]) {
    // create score lookup table as static 3d array
    // calculate Pscore
    for (int i = 0; i < q2; i++) {
        for (int j = 0; j < n2; j++) {
            for (int w = 0; w < n4; w++) {
                table[i][j][w] = get_Pscore(j, w, (double)i + 5);
            }
        }
    }

}

/* A function to compute xcorr.
*/
double xcorr(char *seq, double *frmz, double *frint, int charge, double prec, int n, int *counts, 
		int numMod, Item item, Item item1, mod_ind *mod, results *results1, results *results2, 
		int specNo, char *proteins, int missed, int decoy, int contam, int N_mod, 
		double score_table[][100][140]) {
    int n1 = strlen(seq);
    int n2 = 0; // for counting matched masses
    double res123 = 0;
    double alt_res123 = 0;
    int succss = 0;

    double xcr = 0;
    double xcr_match = 0;

    // 8.11.2022: fr now n_losses will be 0 and defined here. In future could bcast them from main
    int n_losses = 0;

    Item *item_array = NULL; //(Item *)malloc(0);
    int n_item = 0;

    // Update 27.08.2021: need to add information about mod counts
    char probabilities[2048];
    char mod_count[32];

    // 29.10.2022: will implement threshold scoring as in Mascot. Just need the number of candidates
    // for each spectrum
    seq_count++;


    // for sprinf of integers and chars later
    char str22[2048];

    int n4 = 2*n1 - 2;
    int n3 = 0;
    int n5 = 2*n1 - 2; 

    // 30.10.2022: need to refactor here to build separate BST for each isoform
    // This means I need to call BST functions here, not the BST itself
    bst(counts, seq, 42, item, item1, mod, numMod); // 42 is an unnecessary parameter, replaced n, will leave it for now
    
    // debug
    /*if (head->N == 0) {
        printf("BST Head N = %d\n", head->N);
        printf("numMod = %d\n", numMod);
        printf("ptm = %s, %s, %s\n", ptm[0].name, ptm[1].name, ptm[2].name);
    }*/
    
    Item *item_array1 = NULL;
    item_array1 = malloc(head->N*sizeof(Item)); // don't forget to free
    int n_item1 = 0;

    // convert bst to ardered array of Item pointers
    STto_array(head, item_array1, &n_item1);

    int max_match = 0; // this is to store matched 
    int *n2_array = malloc(sizeof(int)*comb_index);
    int *n3_array = malloc(sizeof(int)*comb_index);
    int *n4_array = malloc(sizeof(int)*comb_index);
    if (!n4_array) {
	    printf("Cannot allocate memory for n4_array!\n");
	    exit(1);
    }
    int *n5_array = malloc(sizeof(int)*comb_index);

    for (int i = 0; i < comb_index; i++) {
        n2_array[i] = 0;
        n3_array[i] = 0;
        n4_array[i] = n5;
        n5_array[i] = n5;
    }

    int y_fr[n1];
    int b_fr[n1];

    for (int i = 0; i < n1; i++) {
        y_fr[i] = 0;
        b_fr[i] = 0;
    }
    
    // Now use the tree to search for matches to the spectrum masses. Should account
    // for doubly charged fragments if precursor charge > 2.
    /*for (int i = 0; i < n; i++) {
        xcr = xcr + frint[i];

        succss = interval_search(head, frmz[i] - tol, frmz[i] + tol, &item_array, &n_item, frmz, i, 0, tol, y_fr, b_fr, frint, &n2, &n4, n2_array, numMod, n3_array, n4_array, n5_array);
        //n2 += succss;
	//n3 += succss;
        
        if (charge > 2){
            succss = interval_search(head, (frmz[i]*2 - Pr) - tol, (frmz[i]*2 - Pr) + tol, &item_array, &n_item, frmz, i, 1, tol, y_fr, b_fr, frint, &n2, &n4, n2_array, numMod, n3_array, n4_array, n5_array);
            //n2 += succss;
	    //n3 += succss;
        }


    }*/

    interval_search(head, &item_array, &n_item, frmz,
                        0, tol, y_fr, b_fr, frint, &n2, &n4, n2_array, numMod, n3_array, n4_array,
                        n5_array, n_item1, item_array1, n);
    if (charge > 2) {
        interval_search(head, &item_array, &n_item, frmz,
                        1, tol, y_fr, b_fr, frint, &n2, &n4, n2_array, numMod, n3_array,
                        n4_array, n5_array, n_item1, item_array1, n);
    }

    for (int i = 0; i < n; i++) xcr = xcr + frint[i];

    // 29.09.2022: now calculate intensity coverage
    for (int i = 0; i < n_item; i++) {
        xcr_match += item_array[i]->intensity;
    }

    // 8.11.2022: need to calculate score for each combination because both, n2 and n4 are different
    // and should put this calculation in a separate function at last

    double max_score_n2 = 0;
    double max_score_n3 = 0;
    double tmp_score_n2 = 0;
    double tmp_score_n3 = 0;

    for (int i = 0; i < comb_index; i++) {
        //tmp_score_n2 = get_Pscore(n2_array[i], n4_array[i], q);
        //tmp_score_n3 = get_Pscore(n3_array[i], n5_array[i], q);
        tmp_score_n2 = score_table[(int)(q - 5)][n2_array[i]][n4_array[i]];
        tmp_score_n3 = score_table[(int)(q - 5)][n3_array[i]][n5_array[i]];

        if (tmp_score_n2 > max_score_n2) {
            max_score_n2 = tmp_score_n2;
            n2 = n2_array[i];
            n4 = n4_array[i];
        }
        if (tmp_score_n3 > max_score_n3) {
            max_score_n3 = tmp_score_n3;
            n3 = n3_array[i];
            n5 = n5_array[i];
        }

    }
    
    // debug
    /*printf("n_item1 = %d\n", n_item1);
    printf("comb_index = %d\n", comb_index);
    printf("%s %lf %d %d\n", seq, max_score_n2, n2, n4);
    printf("%s %lf %d %d\n", seq, max_score_n3, n3, n5);*/

    alt_res123 = max_score_n3;
    res123 = max_score_n2;

    // 11.10.2022: only proper y and b matches will be considered for scoring to test if performance better
    //if (alt_res123 > res123 && n_losses != 0) {
    if (n_losses != 0) {
        res123 = alt_res123;
        n2 = n3;
        n4 = n5;
    }

    // Update 16.08.2021: starting refactoring to use new results type
    if (res123 > results1[specNo].Score) {
        results1[specNo].Total = n4;
        results1[specNo].Matched = n2;
	results1[specNo].matched_int = xcr_match/xcr;
        results1[specNo].Score = res123;
        strcpy(results1[specNo].Proteins, proteins);
        strcpy(results1[specNo].Sequence, seq);
        results1[specNo].nMod = N_mod;
        results1[specNo].Missed = missed;
        results1[specNo].Decoy = decoy;
        results1[specNo].Contam = contam;

        if (counts[0] == 0) {
            strcpy(results1[specNo].mod1_probs, "NA");
            strcpy(results1[specNo].mod1, "NA");
        }

        else {
            strcpy(probabilities, "");
            compute_site_prob_open_optimized(probabilities, 0, counts[0], seq, item_array, n_item);
            strcpy(results1[specNo].mod1_probs, probabilities);
            strcpy(results1[specNo].mod1, ptm[0].name);
            sprintf(mod_count, "(%d)", counts[0]);
            strcat(results1[specNo].mod1, mod_count);
            strcpy(mod_count, "");
        }
        if (counts[1] == 0) {
            strcpy(results1[specNo].mod2_probs, "NA");
            strcpy(results1[specNo].mod2, "NA");
        }
        else {
            strcpy(probabilities, "");
            compute_site_prob_open_optimized(probabilities, 1, counts[1], seq, item_array, n_item);
            strcpy(results1[specNo].mod2_probs, probabilities);
            strcpy(results1[specNo].mod2, ptm[1].name);
            sprintf(mod_count, "(%d)", counts[1]);
            strcat(results1[specNo].mod2, mod_count);
            strcpy(mod_count, "");
        }
        if (counts[2] == 0) {
            strcpy(results1[specNo].mod3_probs, "NA");
            strcpy(results1[specNo].mod3, "NA");
        }
        else {
            strcpy(probabilities, "");
            compute_site_prob_open_optimized(probabilities, 2, counts[2], seq, item_array, n_item);
            strcpy(results1[specNo].mod3_probs, probabilities);
            strcpy(results1[specNo].mod3, ptm[2].name);
            sprintf(mod_count, "(%d)", counts[2]);
            strcat(results1[specNo].mod3, mod_count);
            strcpy(mod_count, "");
        }
        strcpy(results1[specNo].Fragments, "");
        strcpy(results1[specNo].Masses, "");
        for (int i = 0; i < n_item; i++) {
            sprintf(str22, "%s,", item_array[i]->fragType);
            strcat(results1[specNo].Fragments, str22);
        }
        for (int i = 0; i < n_item; i++) {
            sprintf(str22, "%lf,", *(item_array[i]->key));
            strcat(results1[specNo].Masses, str22);
        }
	// 11.09.2022: need to free  the items in item_array
        for (int i = 0; i < n_item; i++) {
            free(item_array[i]->key);
	    //free(item_array[i]->comb_ind);
            free(item_array[i]);
        }
        free(item_array);
        STdestroy(head);
	free(z);
	free(n2_array);
	free(n3_array);
	free(n4_array);
	free(n5_array);
	free(item_array1);
        return res123;
    }
    else if (res123 > results2[specNo].Score) {
        results2[specNo].Total = n4;
        results2[specNo].Matched = n2;
	results2[specNo].matched_int = xcr_match/xcr;
        results2[specNo].Score = res123;
        strcpy(results2[specNo].Proteins, proteins);
        strcpy(results2[specNo].Sequence, seq);
        results2[specNo].nMod = N_mod;
        results2[specNo].Missed = missed;
        results2[specNo].Decoy = decoy;
        results2[specNo].Contam = contam;

        if (counts[0] == 0) {
            strcpy(results2[specNo].mod1_probs, "NA");
            strcpy(results2[specNo].mod1, "NA");
        }
        else {
            strcpy(probabilities, "");
            compute_site_prob_open_optimized(probabilities, 0, counts[0], seq, item_array, n_item);
            strcpy(results2[specNo].mod1_probs, probabilities);
            strcpy(results2[specNo].mod1, ptm[0].name);
            sprintf(mod_count, "(%d)", counts[0]);
            strcat(results2[specNo].mod1, mod_count);
            strcpy(mod_count, "");
        }
        if (counts[1] == 0) {
            strcpy(results2[specNo].mod2_probs, "NA");
            strcpy(results2[specNo].mod2, "NA");
        }
        else {
            strcpy(probabilities, "");
            compute_site_prob_open_optimized(probabilities, 1, counts[1], seq, item_array, n_item);
            strcpy(results2[specNo].mod2_probs, probabilities);
            strcpy(results2[specNo].mod2, ptm[1].name);
            sprintf(mod_count, "(%d)", counts[1]);
            strcat(results2[specNo].mod2, mod_count);
            strcpy(mod_count, "");
        }
        if (counts[2] == 0) {
            strcpy(results2[specNo].mod3_probs, "NA");
            strcpy(results2[specNo].mod3, "NA");
        }
        else {
            strcpy(probabilities, "");
            compute_site_prob_open_optimized(probabilities, 2, counts[2], seq, item_array, n_item);
            strcpy(results2[specNo].mod3_probs, probabilities);
            strcpy(results2[specNo].mod3, ptm[2].name);
            sprintf(mod_count, "(%d)", counts[2]);
            strcat(results2[specNo].mod3, mod_count);
            strcpy(mod_count, "");
        }

        strcpy(results2[specNo].Fragments, "");
        strcpy(results2[specNo].Masses, "");
        for (int i = 0; i < n_item; i++) {

            sprintf(str22, "%s,", item_array[i]->fragType);
            strcat(results2[specNo].Fragments, str22);
        }
        for (int i = 0; i < n_item; i++) {
            sprintf(str22, "%lf,", *(item_array[i]->key));
            strcat(results2[specNo].Masses, str22);
        }

        for (int i = 0; i < n_item; i++) {
            free(item_array[i]->key);
	}
	free(item_array);
        STdestroy(head);
	free(z);
	free(n2_array);
	free(n3_array);
        free(n4_array);
        free(n5_array);
	free(item_array1);
        return res123;

    }

    else {
	for (int i = 0; i < n_item; i++) {
	    free(item_array[i]->key);
	    //free(item_array[i]->comb_ind);
	    free(item_array[i]);
	}
        free(item_array);
        STdestroy(head);
	free(z);
	free(n2_array);
	free(n3_array);
        free(n4_array);
        free(n5_array);
	free(item_array1);
        return 0;
    }

}

/* Callback function for sqlite3
   Need to include matched and total ion numbers in the results. Will have to modify the function
   to receive the score and the matched ion number as a pointer to array (18.07.2020).
*/
int callback(sqlite3 *db, char *sql_inst_1, int rowid, double mass, int *counts, int i, char *seq, 
		void *dt, char *proteins, Item item, Item item1, mod_ind *mod, link head, link z, int missed, 
		int decoy, int contam, results *results1, results *results2, int specNo, 
		double score_table[][100][140]) {
    //int i;

    // For debugging
    //for (int ctr = 0; ctr < i; ctr++) printf("%s is %d\n", ptm[ctr].name, ptm[ctr].count);

    spectrum *spec1 = (spectrum *) dt;

    // 02.02.2024: need to copy spectrum because it gets changed for each candidate sequence
    spectrum spec = *(spec1);

    //sqlite3 *db = (sqlite3 *) db1;
    int rc_1;
    double score = 0;
    int N_mod = 0;
    int numMod = 0;
    sqlite3_stmt *ppStmt_1;
    // these will be used to test if string is compatible with modifictions
    int nMod1 = 0;
    int nMod2 = 0;
    int nMod3 = 0;

    /* Update 11.08.2021: starting implementation of open search. Will change the function to accept the database connection
       and delta mass, then search the mod_comb table, populate the ptm array, create mod_ind etc and call xcorr.


       13.08.2021: Need to test if peptide is unmodified: delta mass is within the high-res precursor tol and if so compute score for unmodified       peptide and return
    */

    if (fabs(spec.precMass - mass) < 10*spec.precMass/1000000) {
        score = xcorr(seq, spec.fragMass, spec.fragInt, spec.charge, spec.precMass, spec.actual_size, counts, numMod, item, item1, mod, results1, results2, specNo, proteins, missed, decoy, contam, N_mod, score_table);

        // 29.10.2022: major bug detected. The below code replaces results1 and results2 entries even if score is less than theirs.
	if (score > 0) {
        // it should be if (score > results1[specNo].Score)
	// false alarm, is ok
	//if  (score > results1[specNo].Score) {     
	    results1[specNo].Scan = spec.scan;
            results1[specNo].PrecScan = spec.scanNum;
            results1[specNo].RetTime = spec.retTime;
            results1[specNo].Charge = spec.charge;
            results1[specNo].PrecMass = spec.precMass;
            //return 0;
	//}
	//else if (score > results2[specNo].Score) {
            results2[specNo].Scan = spec.scan;
            results2[specNo].PrecScan = spec.scanNum;
            results2[specNo].RetTime = spec.retTime;
            results2[specNo].Charge = spec.charge;
            results2[specNo].PrecMass = spec.precMass;
            //return 0;
	}
	
            //for (int ctr = 0; ctr < 3; ctr++) fprintf(res, "\tNA");
            //fprintf(res, "\t%d\t%d\t%d\t%d\n", N_mod, missed, decoy, contam);
        
	//else return 0;
        return 0; 
    }
    // Update 31.01.2022: need to use fabs as above? Of course not, there are negative delta masses as well but need to reverse the equation..
   // sprintf(sql_inst_1, "select * from mod_comb where Delta_mass between %lf and %lf;", mass - spec->precMass - 0.02, mass - spec->precMass + 0.02);
    //sprintf(sql_inst_1, "select * from mod_comb where Delta_mass between %lf and %lf;", spec->precMass - mass - 80*(fabs(spec->precMass - mass))/1000000, spec->precMass - mass + 80*(fabs(spec->precMass - mass))/1000000);
    sprintf(sql_inst_1, "select * from mod_comb where Delta_mass between %lf and %lf;", spec.precMass - mass - 5*mass/1000000, spec.precMass - mass + 5*mass/1000000);
    //printf("%s\n", sql_inst_1);
    //exit(1);
    sqlite3_prepare_v2(db, sql_inst_1, -1, &ppStmt_1, 0);
    //exit(0);
    while ((rc_1 = sqlite3_step(ppStmt_1)) != 101) {
        if (rc_1 == 100) {
            // now populate ptm,  test if sequence is compatible with modification and call xcorr
            // for debugging
            //printf("ptm[0].delta_mass as it come from db is %lf\n", sqlite3_column_double(ppStmt_1, 1));
            ptm[0].delta_mass = sqlite3_column_double(ppStmt_1, 1);
            ptm[1].delta_mass = sqlite3_column_double(ppStmt_1, 2);
            ptm[2].delta_mass = sqlite3_column_double(ppStmt_1, 3);
            strcpy(ptm[0].name, (char *)sqlite3_column_text(ppStmt_1, 4));
            strcpy(ptm[1].name, (char *)sqlite3_column_text(ppStmt_1, 5));
            strcpy(ptm[2].name, (char *)sqlite3_column_text(ppStmt_1, 6));
            strcpy(ptm[0].site, (char *)sqlite3_column_text(ppStmt_1, 7));
            strcpy(ptm[1].site, (char *)sqlite3_column_text(ppStmt_1, 8));
            strcpy(ptm[2].site, (char *)sqlite3_column_text(ppStmt_1, 9));
            N_mod = sqlite3_column_int(ppStmt_1, 10);
            ptm[0].count = sqlite3_column_int(ppStmt_1, 11);
            ptm[1].count = sqlite3_column_int(ppStmt_1, 12);
            ptm[2].count = sqlite3_column_int(ppStmt_1, 13);

            // for debugging
            //printf("ptm[0].delta_mass in thread %d is %lf\n", omp_get_thread_num(), ptm[0].delta_mass);

            /* need to check if mod1,mod2,mod3 are different and reshape the mod array accordingly at some point so bst can work properly
               this will be easiest done by changing mod_comb table to contain counts for each mod
               ok, now it's done
            */
            // now check if seq is compatible with the modifications
            /*
               Update 13.08.2021:
               this is not trivial as modifications can have overlapping specificities. Thus we need to account for this by summing up
               somehow. For example:

               let ptm[0].site = "*KN"
                   ptm[1].site = "K"
                   ptm[2].site = "NQ"
               This means that the seq will be compatible only if it has K or N for ptm[0] plus extra K for ptm[1] plus extra N or Q for ptm[2]
               I need a function that will do this test efficiently and this will take some time to develop.
               Algorithm 1:
                 1. Encode the sequence with 0 for each residue.
                 2. For each ptm encode matching sites with corresponding numbers (1, 2 or 3). This will creates 3 strings (or arrays rather).
                 3. "Merge" the 3 arrays into a single final array:
                    for i from 1 to end of string -> if array1[i] is not 0 than final_array[i] = array1[i]; nMod1++;
                    for i from 1 to end of string -> if array2[i] is not 0 and final_array[i] is 0 than final_array[i] = array2[i]; nMod2++;
                    for i from 1 to end of string -> if array3[i] is not 0 and final_array[i] is 0 than final_array[i] = array3[i]; nMod3++;                     This covers the case where specificitites are not overlapping and we simply check as in the main code and proceed.
                    But the check might fail when specificities overlap so we continue
                    if (nMod2 < ptm[1].count) go over the final_aray again and re-encode 0 with 1. At each re-encoding keep track of nMod1,
                    nMod2 and nMod3. test if they became sufficient and if so proceed. If nMod1 or nMod2 become less than necessary
                    undo re-encoding
                    Now do this with ptm[2]....

                    Just realized that I have the algorithm already coded for mod_pep. Should just adapt it and use it here.
                    ... but the mod_pepalgorithm does not allow overlapping specificities so need to change it anyway

                    1. Step through the string and check if residue is in the specificity string of the 3 ptms (2 nested loops:
                       outern is over ptms and inner is over the string)
                       increment a counter variable. This is the length of the multiset. From the example specificities above
                       let seq is *LGGKNPPQQR and let counter is y then we have in the two loops:
                       for ptm[0]:  y becomes 3
                       for ptm[1]:  y becomes 4
                       for ptm[2]:  y becomes 7 so the multiset has 7 elements

                       now we step through the sequence first and through ptms in an inner loop and populate the multiset array:
                       0010222
                       Now sort it in lexicogrpahical order: 0001222 so the peptide has 3 potential ptm[0] sites, 1 ptm[1] site
                       and 3 ptm[2] sites but there is overlap and how should we deal with this? If the delta mass suggest 3 ptm[0]
                       modifications this means no ptm[1] is possible. Should change the algorithm to correct multiset for overlaps
                       but how?

                       Let's explore sets and algebra: if the set of ptm[i] sites and ptm[j] sites are overlaping then we neeed the
                       union to form multiset but the union. What the above does is not the union (it only generates the union for
                       non-overlapping specificities). It doubles the intersection so it is the union + intersection. Thus the proper
                       algorithm could either directly find the union, or, if it is more efficient, use the above, find the intersection
                       and subtract it from the multiset.


                       So. let's try a new algrithm:
                       Step over the seq and ecode three arrays for ptm1,2, and 3 (0 non-site, ptm[0] = 1 etc):
                       for ptm[0]: 10001100000 nMod1 = 3
                       for ptm[1]: 00002000000 nMod2 = 1
                       for ptm[2]: 00000300330 nMod3 = 3

                       now we make a union and encode the element that has overlaps with the number that has least counts
                       and correct the nMod accordingly:
                       union is: 10002(nMod1 becomes 2)1(nMod3 becomes2)00330 so this sequence can accomodate the up to 1 ptm[1] and
                       either 2 ptm[2] and 3 ptm[3] or 3 ptm[2] and 2 ptm[3] but how to implement the "either" part

                       the union is: 1000(1 or 2)(1 or 3)0030 we should then generate combinations perhaps and test if any combination
                       is compatible with the delta mass the combinations are:
                           {1,1,1,3,3} no ptm2
                           {1,1,2,3,3}
                           {1,1,3,3,3} no ptm2
                       so the algorithm should generate the combinations and test if there is at least one that is compatible with the
                       ptm array.

                       so the multiset should have 5 elements and we generate combinations and test if any satisfies the conditions:
                       we can modify the mod_pep multiset algorithm:
                       First we determine the smalest nMod and this is nMod2 = 1. If corresponding ptm[1].counter is not 0 we need to
                       take care for this first so we assign at least ptm[1].counter elements of the multiset to this ptm:
                       multiset = {2....} then we rucurr and find the next smallest and assign accordingly and so on until we fill the multiset

                       actually we should start with the smallest ptm[].count, that maks more sense then assign necessary elements if the
                       nMod allows, if not return without processing. If it does then proceed with the next in order

                       Update 30.01.2022: still not working with overlapping site specificities. Should give it couple of days and if not
                       fixed should write up and submit only scorer_omp and leave combinatorial search for next paper. Let's see...

                       could if else if work? Lets try to make a function:

                       it will take ptm array and sequence and return 1 if they are compatible and 0 if not.
                       If any two ptms have overlapping specificities what should the function do? Concatenate site strings? 
                       no, let's try something different:
                       
                       declare 3 new variable: overlap_12, overlap_13, and overlap_23. Initialise them to 0

                       1. Loop through seq and for each char in seq:
                          test if seq[w] matches ptm[0].site
                            If yes, test if it matches ptm[1].site
                              if not, test if it matches ptm[2].site
                                if no increment nMod1
                                if yes increment overlap_13 variable
                              if yes increment overlap_12
                            If not, test if it matches ptm[2].site
                              if yes test if it matches ptm[3].site
                                if not increment nMod2
                                if yes increment overlap_23
                              if not test if it matches ptm[3].site
                                if yes increment nMod3
                       2. At the end of the loop we have nMod and overlap variables with counts of available sites. 
                          we need to somehow use them to determine if seq and ptm are compatible: for example, 
                          for ptm[0] we will have total available sites nMod1 + overlap_12 + overlap_13
                          for ptm[1] we will have total available sites nMod2 + overlap_12 + overlap_23
                          for ptm[2] we will have total available sites nMod3 + overlap_13 + overlap_23

                          Now we need to test all possibilities:
                       
                          if nMod1 >= ptm[0].count then seq is compatible with ptm[0]
                            if not we test if nMod1 + overlap_12 + overlap_13 >= ptm[0].count
                              if yes then ptm[0] is compatible with seq but on condition that ptm[1] and ptm[2] do not need the overlaps
                                therefore we now test if nMod2 >= ptm[1].count
                                  if yes then ptm[1] and seq are compatible
                                    then we test if nMod3 >= ptm[2].count
                                      if yes ptm and seq are compatible and the function terminates with return value of 1
                                      if not we test if nMod3 + overlap_13 + overlap_23 >= ptm[2].count
                                        if not the function terminates and returns 0
                                        if yes ptm[0] and ptm[2] need to be considered together: we need to use set theory for this...
                                          will resume tomorrow 31.01.202
     
                                          after much thinking decided that for now it will only look for non-overlapping specificities
                                          will modify mod_DB1.c to make sure it does (today is 31.01.2022)
        
                                  if not we test if nMod2 + overlap_12 + overlap_23 >= ptm[1].count
                                    Now if yes we need to consider ptm[0] and ptm[1] together
                                      
                                    if not the function terminates and returns 0 
                              if not the function terminates and returns 0
                              
            */
            for (int w = 0; w < strlen(seq); w++) {
                if (strchr(ptm[0].site, seq[w])) nMod1++;
                if (strchr(ptm[1].site, seq[w])) nMod2++;
                if (strchr(ptm[2].site, seq[w])) nMod3++;
            }
            if (nMod1 >= ptm[0].count && nMod2 >= ptm[1].count && nMod3 >= ptm[2].count && strlen(seq) > 1) {
                for (int j = 0; j < i; j++) {
                    if (ptm[j].count > 0) numMod++;
                    counts[j] = ptm[j].count;
                }
                score = xcorr(seq, spec.fragMass, spec.fragInt, spec.charge, spec.precMass, 
				spec.actual_size, counts, numMod, item, item1, mod, results1, results2, 
				specNo, proteins, missed, decoy, contam, N_mod, score_table);

                // 29.10.2022: major bug, see above—no it is not
		if (score > 0) {
		//if  (score > results1[specNo].Score) {
                    //for (int ctr = 0; ctr < 3; ctr++) fprintf(res, "\t%s", ptm[ctr].name);
                    //fprintf(res, "\t%d\t%d\t%d\t%d\n", N_mod, missed, decoy, contam);
                    results1[specNo].Scan = spec.scan;
                    results1[specNo].PrecScan = spec.scanNum;
                    results1[specNo].RetTime = spec.retTime;
                    results1[specNo].Charge = spec.charge;
                    results1[specNo].PrecMass = spec.precMass;
                    
		//}
		//else if  (score > results2[specNo].Score) {
                    results2[specNo].Scan = spec.scan;
                    results2[specNo].PrecScan = spec.scanNum;
                    results2[specNo].RetTime = spec.retTime;
                    results2[specNo].Charge = spec.charge;
                    results2[specNo].PrecMass = spec.precMass;
		    
                }
            }
                
		// now need to reset the ptm array for the next candidate
                // update 15.08.2021: counts[] also need resetting for the unmodified peptide search to work
                for (int j = 0; j < 3; j++) {
                    ptm[j].count = 0;
                    counts[j] = 0;
                    // 02.02.2024: also name and site should be reset
                    strcpy(ptm[j].name, "");
                    strcpy(ptm[j].site, "");

                }
                numMod = 0;
                nMod1 = 0;
                nMod2 = 0;
                nMod3 = 0;
		N_mod = 0;
                
            
        }
    }
    sqlite3_finalize(ppStmt_1);

    /*
    if (strlen(seq) > 1) {
        //printf("Sequence is %s.\n", argv[0]);
        score = xcorr(seq, spec->fragMass, spec->fragInt, spec->charge, spec->precMass, spec->actual_size, counts, i, item, item1, mod);
        // For debugging
        //printf("Score is %lf\n", score);

        if (score > 0) {
            fprintf(res, "\b \t%d\t%d\t%lf\t%s\t%d\t%lf\t%lf\t%s", spec->scan, spec->scanNum, spec->retTime, seq, spec->charge, spec->precMass, score, proteins);
            for (int ctr = 0; ctr < i; ctr++) fprintf(res, "\t%d", counts[ctr]);
            fprintf(res, "\t%d\t%d\t%d\t%d\n", nMod, missed, decoy, contam);
        }
    }*/

    // Update 01.08.2020: Decoy analysis. Will generate a reverse sequence and analyse it as above.
    // not working so going back to original. Decoy sequences are now in db
    /*char revSeq[strlen(seq)];
    char s[strlen(seq)];
    strcpy(revSeq, seq);
    revSeq[strlen(seq)-1] = '\0';
    strcpy(s, strrev(revSeq));
    strcat(s, seq + strlen(seq) - 1); // this adds the C-terminal K or R

    // Now proceed with analysis of reverse sequence
    score = 0;

    char decoy[2048] = "decoy_:";

    strcat(decoy, proteins);
    score = xcorr(s, spec->fragMass, spec->fragInt, spec->charge, spec->precMass, sizeof(spec->fragMass)/sizeof(double), counts, i, item, item1, mod);

    if (score > 0) {
        fprintf(res, "%d\t%d\t%lf\t%s\t%d\t%lf\t%lf\t%s", spec->scan, spec->scanNum, spec->retTime, s, spec->charge, spec->precMass, score, decoy);
        for (int ctr = 0; ctr < i; ctr++) fprintf(res, "\t%d", counts[ctr]);
        fprintf(res, "\n");
    }

 */

    return 0;
}

int scorer(spectrum spec, results *results1, results *results2, int specNo, sqlite3 *db1, Item item, 
		Item item1, PTM *ptm, double precTol, double tol, double score_table[][100][140]) {

    // update 29.08.2021: need to declare and assign precTol and tol here or broadcast from master process
    // update 22.09.2021: now changing this to come from master in a broadcast
    //double precTol = 500;
    //double tol = 0.02;

    mod_ind *mod = NULL; // This works, definition with malloc does not.

    // Array to hold modification str.
    //char modStr[1024];
    char seq[1024];
    char proteins[4096]; // Had to increase it to 4096 as it was giving error running on mac os.
    //char *token;
    int i = 3; // now fixed to 3


    // Create counts[i] array to hold number of modifications of each type..
    // this should be passed to callback and populated there
    int counts[i];
    for (int w = 0; w < i; w++) counts[w] = 0;

    // Open sqlite database.

    sqlite3_stmt *ppStmt;
    //sqlite3_stmt *ppStmt_2 = NULL;

    //char *zErrMsg = 0; // Not used.
    int rc;

    double lower;
    double upper;
    int ct = 0;

    char sql_inst[1024];
    char sql_inst_1[1024];

    //Item item = NULL;
    //Item item1 = NULL;

    // now changing the while loop to parallel for



            // Create item pointers and allocate memory.
            /*item = malloc(sizeof(*item));  // This might need to go to scorer_mod.c and be passed as argument.
            if (!item) {
                printf("Error allocating memory for item!\n");
                exit(1);
            }

            // for debugging
            //printf("Created item in thread %d address is %x\n", omp_get_thread_num(), &item);

            item1 = malloc(sizeof(*item1)); // This might need to go to scorer_mod.c and be passed as argument.
            if (!item1) {
                printf("Error allocating memory for item!\n");
                exit(1);
            }

            // for debugging
            //printf("Created item1 in thread %d address is %x\n", omp_get_thread_num(), &item1);

            ptm = (PTM *) malloc(3 * sizeof(PTM));
            if (!ptm) {
                printf("Error allocating memory for ptm array!\n");
                exit(1);
            }*/

            //printf("Created ptm in thread %d address is %x\n", omp_get_thread_num()), &ptm;

            //fseek(mms_ptr, sizeof(spectrum) * specNo, SEEK_SET);
            //fread(&spec1, sizeof(spectrum), 1, mms_ptr);

            lower = spec.precMass - precTol;
            upper = spec.precMass + precTol;

            sprintf(sql_inst, "select * from peptides where mass between %lf and %lf;", lower, upper);
            //printf("Searching for %lf...\n", spec.precMass); // For troubleshooting only
            // Should use prepare/step instead of exec.
            // rc = sqlite3_exec(db, sql_inst, callback, &spec1, &zErrMsg);

            rc = sqlite3_prepare_v2(db1, sql_inst, -1, &ppStmt, 0);

            while ((rc = sqlite3_step(ppStmt)) != 101) {
                if (rc == 100) {

                    /* Read the sequence and modifications and call scoring function.
                       The columns that we need are:
                       0. Rowid (pept id), col 0;
                       1. Mass, col 1;
                       2. 2 to 2 + i - 1. Mod counts 0 to i - 1 (the number of elements in ptm array);
                       3. 2 + i. Sequence;

                    */
                    int rowid = sqlite3_column_int(ppStmt, 0);
                    double mass = sqlite3_column_double(ppStmt, 1); // This will be needed later when we use delat mass.
                    //int counts[i]; // This should only be created once outside the loop.
                    //for (int ctr = 0; ctr < i; ctr++) counts[ctr] = sqlite3_column_int(ppStmt, ctr + 2);

                    // This made me spend a day with valgrind because ptm[].count were not updated.
                    // Do I need counts—I can read the info from ptm.
                    //for (int ctr = 0; ctr < i; ctr++) ptm[ctr].count = sqlite3_column_int(ppStmt, ctr + 2);

                    strcpy(seq, (char *) sqlite3_column_text(ppStmt, 0));

                    // For debuging only
                    //printf("Sequence is %s\n", seq);

                    // Need to parse the string to numbers later when summarising results.
                    strcpy(proteins, (char *) sqlite3_column_text(ppStmt, 2));
                    int missed = sqlite3_column_int(ppStmt, 3);
                    int decoy = sqlite3_column_int(ppStmt, 4);
                    int contam = sqlite3_column_int(ppStmt, 5);

                    callback(db1, sql_inst_1, rowid, mass, counts, i, seq, &spec, proteins, item, item1, mod, head, z,
                             missed, decoy, contam, results1, results2, specNo, score_table);
                    // head and z are extern so they don't need to be passed to call back and xcorr. should fix this at some point

                    /*
                    int count_sum = 0;
                    for (int ctr = 0; ctr < i; ctr++) {
                        count_sum += counts[ctr];
                        //printf("%s is %d\n", ptm[ctr].name, ptm[ctr].count);
                    }
                    if (count_sum <= maxMod) { // Need to limit total number of modifications to be considered to avoid hanging the program.
                        callback(rowid, counts, i, seq, &spec1, proteins, item, item1, mod, head, z, missed, decoy, contam, count_sum); // deleted ptm for debugging
                    }*/
                    //else printf("Peptide is not modified.\n"); // For debugging only
                    //callback(rowid, counts, i, seq, &spec1, proteins, item, item1, mod, head, z); // deleted ptm for debugging
                } else {
                    printf("Error executing sql statement: %d.\n", rc);
                    exit(1);
                }
            }

            ct++;
            //ct1 += ct;
            sqlite3_finalize(ppStmt);




            //free(ptm);
            //free(item);
            //free(item1);

            /*
             * Update 17.09.2021: will move declaration of item, item1
             * and ptm to main and make them automatic variables to avoid
             * malloc and free
             */








    //sqlite3_finalize(ppStmt);
    //sqlite3_close(db);

    //free(ptm);
    //free(item);
    //free(item1);
    //free(results1);
    //free(results2);

    //free(mod);
    //free(z);
    // Test with madeup data.
    //double frmz[] = {499.2,823.2,1272.4,331.16, 234.1, 216.1, 762.37, 1117.5};
    //double frint[] = {2500.0,3250.0, 3250.0, 10275.0, 315.0, 1225.0, 4500.0, 2000.0};
    //double prec = pepMass(argv[1]);
    //int n_fr = sizeof(frmz)/sizeof(frmz[0]);
    //double score = xcorr(argv[1], frmz, frint, prec, n_fr);
    //printf("Score for %s is %lf.\n", argv[1], score);

    mpfr_free_cache(); // To prevent memory leaks

    // for debugging
    //printf("%s\n", results1[specNo].Sequence);

    return 0;
}
