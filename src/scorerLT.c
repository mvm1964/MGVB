/*
  Minimal MS/MS search engine in C 
  Usage: scorer_2020 <mms file> <db file> <config file> <maxMod>
  Metodi V. Metodiev, 2019, UK     

  Compilation:
  gcc -g -Wall scorer_mod_2020.c -o scorer_mod_2020_linux -lgmp -lmpfr -lsqlite3 pepMass.c Item.c bst.c

  Update 12 July 2019: starting impementation of modified peptide search.

  Algorithm notes:
  1. Uses pointers to mod_ind and ptm arrays as in comb_kN_new.c.
  2. Retrieves modifications info, peptide and protein IDs and sequence from pep_mod table.
  3. For each predicted fragment ion generate acceptable combinations of
     masses derived from modified sequences. Acceptable means all combinations that satisfy
     the limits set by counts[]. Cognate y and b ions are considered in this proces. For 
     example a doubly phosphorylated peptide with a sequence GESPQSPSPR cannot have a unmodified
     y5 and beyond ions and also cannot have unmodified b ions beyond b6. On the other hand, the 
     peptide can produce singly and doubly phosphorylated y5 fragments and so on. 
     This could be done by calling a function similar to comb_kN_new.c with fragment sequence and
     mod info, It will generate modified fragments masses list, which can be used for scoring.
     For exampe:
         Let sequence of peptide is GESPQSPSPR and it is doubly phosphorylated.
         Consider the SPSPR fragment (which is y5). Acceptable modifications are:
            1 and 2 phosphorylations—0 is not acceptable as there is only one modifiable
            residue in the cognate b6 peptide..    
  4. Computes score based on ALL possible modified fragment masses INCLUDED.
     Thus it covers mixturs of peptide precursor ions modified on alternative sites.
  5. For a selected set of hits generates all combinations of uniquely modified sequences
     and computes corresponding individual scores.
  6. Uses the individual scores to calculate localisation probabilities. 
     
  Update 01.08.2020: starting the implementation of decoy searches within the scorer program. Instead
  of puting decoy sequnces in the sqlite tabes, scorer will generate the reverse sequences on the fly 
  and seqrch against the MS/MS peaks.
  Update 03.08.2020: On the fly analysis of reverse decoy peptide sequences generates high scoring sequences
  and does not work for filtering. Possibly this is only vaid for modified peptides but still, need to switch 
  bck to using reversed whole-protein sequences and implement this in the digestor program.  

  Update 15.08.2020: starting the parser function that should read from the config file and parse

  - precTol
  - tol
  - FASTA files
  - maxMod
  - Variable modifications
  - Fixed modifications

  Will use fgets() as in parseMS
   
  28.09.2022: will implement intensity scoring and phospho neutral loss to improve fdr
  will go like this: sum up the intensity of all peaks matched to theoretical ions. Will record the 
  intensities in the item_array, items should have a field intensity. 
  No, simples is to define a variable, success and assign it a value of 1 interval_search is
  successful or 0 if not. Then add frint[i]*success to a variabe xcr_match.

  08.10.2022: major refactoring. Instead of results1 and 2 will pass only a single results struct to xcor.
  Then in call back will use it to populate results1 and results2

  8.11.2022: starting major refactoring:
  n2_array will hold the matched fragments for each combination with no neutral losses considered
  n3_array will hold the matched fragments with neutral losses considered
  n4_array will hold total expected fragments with no neutral losses considered
  n5_array will hold total with neutral losses considered.
*/

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
#include "compute_site_prob_optimized.h"
#include <omp.h>
#include "parse_spectrum.h"
#include "arraySearch2023.h"
#include "mass_recalibration.c"

/* update 02.09.2021: will optimize performance by 
   1. reduce writes to file
   2. move db to memory
   for 1 will make results structs similar to scorer_open_mpi but simpler 
   
   03.09.2021: after much work it does not speed up when only top and next best are saved to file.
   Will now try moving db to memory
*/

/*
typedef struct Results {
    double score;
    char data[4096];
} results;
*/

// This function reverses strings (copied it from the web).
/*char *strrev(char *str) {
      char *p1, *p2;

      if (! str || ! *str)
            return str;
      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}*/


// calculates p-score
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

    // for debuging
    //printf("n2, n4, q1, sore: %d, %d, %d, %lf\n", n2, n4, q1, res123);
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

// Parser function for parsing run parameters
int parse_params(char *config_file, double *prTol, double *frTol, char fasta[5][2048], 
		int *maxMod, double *q_start, double *q_max, int *n_losses, 
		int *check_alt_iso, double *prTol_1st, double *min_score, int *max_data) {

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
            sscanf(str, "_precTol\t%lf", prTol);
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
	else if (strstr(str, "_qmax")) {
            sscanf(str, "_qmax\t%lf", q_max);
        }
        else if (strstr(str, "_qstart")) {
            sscanf(str, "_qstart\t%lf", q_start);
        }
        else if (strstr(str, "_nLosses")) {
            sscanf(str, "_nLosses\t%d", n_losses);
        }
	else if (strstr(str, "_check_alt_iso")) {
            sscanf(str, "_check_alt_iso\t%d", check_alt_iso);
        }
	else if (strstr(str, "_precTol_1st")) {
            sscanf(str, "_precTol_1st\t%lf", prTol_1st);
        }
	else if (strstr(str, "_min_score")) {
            sscanf(str, "_min_score\t%lf", min_score);
	}
	else if (strstr(str, "_max_data")) {
            sscanf(str, "_max_data\t%d", max_data);
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
double precTol = 10;
double tol = 0.6;
double q_max = 7.0;
double q_start = 1.0;
int n_losses = 0;
int check_alt_iso = 0;
double precTol_1st = 20;
double min_score = 120;
int max_data = 500;

/* A function to compute xcorr.
   
   Now it considers modifications (see update above). It should do something like:

   If there are modification:

   1. For each fragMass: look in the sequence and mod info and decide what masses
      to compare to the spectral data. Should generate combinations, test
      if corresponding masses match spectral peaks and update the counter.
      Combinations: 
      - First, for each type of modification determine the total number of modifiable 
        residues (N[], should be an array of size i).
      - For a specific y fragment determine the number of modifiable residues of given 
        type (n1[j]).
      - Let the required number for a modification be n (from the data returned from mod-pep query).
      - Now:
            if (n1[j] = 0) no modifications for this fragment;
            else if (N[j] - n1[j] < n) combine from n - (N[j] - n1)[j] to n1[j];
            else combine from 0 to n1[j];
            cognate fragment masses are computed as usual: 
                pepMass - fragMass + Pr + Pr;  
   2. Keep account of what peptide is top scorer.

   If top scorer is modified peptide:

   3. Estimate localisation probabilities for the top scorer.
  
   If top scorer is unmodified proceed as previously.

   Update 09.09.2019

      It might be simpler to go to individual combinations as geenrated by comb_kN.c and
      then sum up the matches to allow for isoforms. Given that com_kN is already implemented
      this will speed up development and will also allow easier computation of individual scores for
      isoforms.

      to do this we will need code that uses  
   
   Update May 2020: now we use a bst function which generates a search tree for all possible combinations
   of modifications positions and corresponding y and b fragments. scorer_mod now needs to use this tree to
   check each spectrum peak for matching and output the corresponding data to a file. It also needs to use the
   Bayesian network function to compute site probabilities.


   Update 02.08.2021: need to refactor xcorr and calback to be able to only write results rows with score > of certain limit

*/
double xcorr(char *seq, double *frmz, double *frint, int charge, double prec, int n, 
		int *counts, int numMod, Item item, Item item1, mod_ind *mod, 
		results *xcr_results, char *str3, double q1, int n_item1, Item *item_array1, double score_table[][100][140]) {
    int n1 = strlen(seq);
    int n2 = 0; // for counting matched masses
    int n3 = 0; // for counting only y and b fragments 06.10.2022
    int n5 = 2*n1 - 2; // for alternative scoring 06.10.2022
    double xcr = 0;
    double xcr_match = 0;
    double res123 = 0;
    double alt_res123 = 0;
    
    /* Update 15.07.2020: create item_array and n_item vars. Should not forget to free(item_array) 
       because interval_search() realloc memory to it.
    */
    Item *item_array = NULL; //(Item *)malloc(0);
    int n_item = 0;


    // Update 19.07.2021: p should be a function of fragment mass tolerance and the number of peaks in the spectrum
    
    // bst should now be created in callback
    //bst(counts, seq, 0, item, item1, mod, numMod);
    //STprint(head);
    //exit(1);
    
    int *n2_array = malloc(sizeof(int)*comb_index);
    int *n3_array = malloc(sizeof(int)*comb_index);
    int *n4_array = malloc(sizeof(int)*comb_index);
    int *n5_array = malloc(sizeof(int)*comb_index);
    for (int i = 0; i < comb_index; i++) {
        n2_array[i] = 0;
	n3_array[i] = 0;
        n4_array[i] = n5;
	n5_array[i] = n5;
    }

    int n4 = 2*n1 - 2;
    int y_fr[n1];
    int b_fr[n1];
    for (int i = 0; i < n1; i++) {
        y_fr[i] = 0;
	b_fr[i] = 0;
    }

    // Now use the tree to search for matches to the spectrum masses. Should account 
    // for doubly charged fragments if precursor charge > 2. 
   /* for (int i = 0; i < n; i++) {
        //if (frmz[i] == 0 || frint[i] == 0) continue;

        xcr = xcr + frint[i];

	interval_search(head, frmz[i] - tol, frmz[i] + tol, &item_array, &n_item, frmz, 
			i, 0, tol, y_fr, b_fr, frint, &n2, &n4, n2_array, numMod, n3_array, n4_array, n5_array);
        if (charge > 2){
            interval_search(head, (frmz[i]*2 - Pr) - tol, (frmz[i]*2 - Pr) + tol, &item_array, &n_item, frmz, i, 
			    1, tol, y_fr, b_fr, frint, &n2, &n4, n2_array, numMod, n3_array, n4_array, n5_array);
        }
    }*/

    interval_search(head, &item_array, &n_item, frmz,
                        0, tol, y_fr, b_fr, frint, &n2, &n4, n2_array, numMod, n3_array, n4_array,
                        n5_array, n_item1, item_array1, n);
    if (charge > 2){
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
        //tmp_score_n2 = get_Pscore(n2_array[i], n4_array[i], q1);
        //tmp_score_n3 = get_Pscore(n3_array[i], n5_array[i], q1);
        tmp_score_n2 = score_table[(int)(q1 - 5)][n2_array[i]][n4_array[i]];
        tmp_score_n3 = score_table[(int)(q1 - 5)][n3_array[i]][n5_array[i]];

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
    
    alt_res123 = max_score_n3;
    res123 = max_score_n2;

    // 11.10.2022: only proper y and b matches will be considered for scoring to test if performance better
    //if (alt_res123 > res123 && n_losses != 0) {
    if (n_losses != 0) {
        res123 = alt_res123;
	n2 = n3;
	n4 = n5;
    }
    
    // now populate xcr_results
    if (res123 > xcr_results->score) {
	xcr_results->score = res123;
        if (n2 > 0) {
            strcpy(xcr_results->data, "");
            sprintf(str3, "%d\t%d\t%lf\t%lf\t", n4, n2, xcr_match/xcr, q1);
            strcat(xcr_results->data, str3);
            for (int ctr = 0; ctr < numMod; ctr++) {
                if (counts[ctr] == 0) strcat(xcr_results->data, "NA\t");
                else compute_site_prob_optimized(xcr_results, str3, ctr, counts[ctr], seq, item_array, n_item); // need to modify compute_site_prob to strcat
            }
        }
        for (int i = 0; i < n_item; i++) {
            sprintf(str3,"%s,", item_array[i]->fragType);
            strcat(xcr_results->data, str3);
        }
        //fprintf(res, "\t");
        for (int i = 0; i < n_item; i++) {
            sprintf(str3,"%lf,", *(item_array[i]->key));
            strcat(xcr_results->data, str3);
        }

        // free the item_array
        for (int i = 0; i < n_item; i++) {
            free(item_array[i]->key);
            free(item_array[i]);
        }
        free(item_array);
        //STdestroy(head);
	free(n2_array);
	free(n3_array);
	free(n4_array);
	free(n5_array);
        return res123;

    }
    else {
        // free the item_array
        for (int i = 0; i < n_item; i++) {
            free(item_array[i]->key);
            free(item_array[i]);
        }
        free(item_array);
        //STdestroy(head);
	free(n2_array);
	free(n3_array);
        free(n4_array);
	free(n5_array);
        return res123;
    }
}

/* Callback function for sqlite3
   Need to include matched and total ion numbers in the results. Will have to modify the function 
   to receive the score and the matched ion number as a pointer to array (18.07.2020).
*/
int callback(int rowid, int *counts, int i, char *seq, void *dt, char *proteins, Item item, 
		Item item1, mod_ind *mod, int missed, int decoy, int contam, int nMod, 
		results *results1, results *results2, char *str3, double score_table[][100][140]) {
    //int i;

    // For debugging
    //for (int ctr = 0; ctr < i; ctr++) printf("%s is %d\n", ptm[ctr].name, ptm[ctr].count);

    //double q_max = q;
    spectrum *spec = (spectrum *) dt;
    results xcr_results;
    xcr_results.score = 0;
    strcpy(xcr_results.data, "");

    // update 07.10.2022: will start the implementation of varying q: for low res q 5, 6, and 7 will be tested
    // for high res q 8, 9, 10 will be tested
    // to do this will have to use the interval_fill and other functions from parseMS.c
    
    // Update 8.12.2023: major flaw detected: should not build BST for
    // every spectrum in this loop. Should first build BST then create
    // different spectra and pass them along with the BST to a refactored
    // xcorr!!!! This should speed up greately!
    // Why am I passing spec->actual_size to bst????
    //bst(counts, seq, 0, item, item1, mod, i);    

    // debug
    /*printf("Created BST!\n");
    printf("Sequence is %s\n", seq);
    STprint(head);
    STdestroy(head);
    exit(1);*/

    Item *item_array1 = NULL;
    item_array1 = malloc(head->N*sizeof(Item)); // don't forget to free
    int n_item1 = 0;

    // convert bst to ardered array of Item pointers
    STto_array(head, item_array1, &n_item1);

    // 06.01.2024: compute mass error
    double ptm_mass = 0;
    for (int j = 0; j < i; j++) {
        ptm_mass += counts[j]*ptm[j].delta_mass;
    }

    double t_mass = pepMass(seq) + ptm_mass + Pr;
    double MassError = 1000000*(spec->precMass - t_mass)/(t_mass);

    double max_score = 0;
    for (int q1 = q_start; q1 <= q_max; q1++) {
        // continue tomorrow: 08.10.2022
        spectrum spec1 = *(spec);
	for (int t = 0; t < spec1.actual_size; t++) {
            spec1.fragMass[t] = 0;
	    spec1.fragInt[t] = 0;
	}

        int wi[21];
        //for (int j = 1; j <= 21; j++) wi[j] = 0;
	for (int j = 0; j < 21; j++) wi[j] = 0;
        //int q1 = atoi(argv[1]);
        double mass;
        double intensity;
        int interval = 0;

        // now use the spec to create spectra with various q
	// Update 8.12.2023: major flaw detected: should not build BST for 
	// every spectrum in this loop. Should first build BST then create 
	// different spectra and pass them along with the BST to a refactored
	// xcorr!!!! This should speed up greately!
        for (int t = 0; t < spec->actual_size; t++) {
            mass = spec->fragMass[t];
            intensity = spec->fragInt[t];
            // Fill intervals here.
            if (mass > 2000) interval = 21;
            else interval = ceil(mass/100);
            intervalFill(mass, intensity, &spec1, wi, interval, q1);
            wi[interval]++;
        }
        sort_masses(&spec1);
        remove_0(&spec1);
        double score = 0;
        if (strlen(seq) > 1) {
        //printf("Sequence is %s.\n", argv[0]);
            score = xcorr(seq, spec1.fragMass, spec1.fragInt, spec1.charge, spec1.precMass, 
			    spec1.actual_size, counts, i, item, item1, mod, &xcr_results, 
			    str3, q1, n_item1, item_array1, score_table);
            if (score > max_score) {
                max_score = score;
	    }
        }
    }
    
    if (max_score > results1->score) {
        results1->score = max_score;

	strcpy(results1->data, xcr_results.data);
        sprintf(str3, "\b \t%d\t%d\t%lf\t%s\t%d\t%lf\t%lf\t%lf\t%s", spec->scan, spec->scanNum, 
			spec->retTime, seq, spec->charge, spec->precMass, MassError, max_score, proteins);
        strcat(results1->data, str3);
        for (int ctr = 0; ctr < i; ctr++) {
            sprintf(str3, "\t%d", counts[ctr]);
            strcat(results1->data, str3);
        }
        sprintf(str3, "\t%d\t%d\t%d\t%d\t", nMod, missed, decoy, contam);
        strcat(results1->data, str3);
    }
    else if (max_score > results2->score) {
        results2->score = max_score;
	strcpy(results2->data, xcr_results.data);
        sprintf(str3, "\b \t%d\t%d\t%lf\t%s\t%d\t%lf\t%lf\t%lf\t%s", spec->scan, spec->scanNum, 
			spec->retTime, seq, spec->charge, spec->precMass, MassError, max_score, proteins);
        strcat(results2->data, str3);
        for (int ctr = 0; ctr < i; ctr++) {
            sprintf(str3, "\t%d", counts[ctr]);
            strcat(results2->data, str3);
        }
        sprintf(str3, "\t%d\t%d\t%d\t%d\t", nMod, missed, decoy, contam);
        strcat(results2->data, str3);
    }

    //STdestroy(head);
    free(item_array1);
    return 0;
}

int main(int argc, char **argv){
   
    // testing the score lookup table
    static double score_table[6][100][140];
    make_score_table(100, 140, 6, score_table);
    printf("Testing the score lookup table: %lf, %lf, %lf\n",  score_table[0][20][50],  
		    score_table[3][40][42], score_table[5][20][22]);
    //exit(1);


    // debug
    /*double test_score = 0;
    test_score = get_Pscore(atoi(argv[1]), atoi(argv[2]), atof(argv[3]));
    printf("Test score = %lf\n", test_score);
    exit(1);*/

    char fasta[5][2048]; // up to 5 fasta files can be used
    maxMod = 0;

    int first_search = atoi(argv[4]); // should call scorer with 1 for first search and 0 for main search in argv[4]

    char resFile[128] = {};
    strcat(resFile, argv[1]);
    if (first_search == 1) strcat(resFile, ".1st");
    strcat(resFile, ".txt");

    // declare and initialize results1 and results2
    results results1;
    results results2;
    results1.score = 0;
    strcpy(results1.data, "");
    results2.score = 0;
    strcpy(results2.data, "");
    char str3[4096];
    strcpy(str3, "");

    //int N_spectra = atoi(argv[4]);
    int N_spectra = 0;
    int res_count = 0; // 9.12.2023: to count results records
    char str22[1024]; // for forming reslts strings

    // this is to keep track of candidate sequences for the mascot-like scoring
    // has to be private
    int seq_count = 0;

    // Parse config file
    int fasta_cntr = parse_params(argv[3], &precTol, &tol, fasta, &maxMod, &q_start, &q_max, 
		    &n_losses, &check_alt_iso, &precTol_1st, &min_score, &max_data);

    Item item = NULL;
    Item item1 = NULL;
    mod_ind *mod = NULL; // This works, definition with malloc does not.

    // File pointer to read config info.
    FILE *cnfg;
    if ((cnfg = fopen(argv[3], "r")) == NULL) {
        printf("Error opening config file for reading.\n");
        exit(1);
    }

    // Array to hold modification str.
    char modStr[1024];
    char seq[1024];
    char proteins[4096]; // Had to increase it to 4096 as it was giving error running on mac os.
    int i = 0;

    // update 05.09.2021: all dynamically allocated arrays that are to be private should be converted to static
    // omp does not know how much memory to reserve for them and treats them as shared
    
    PTM *ptm_tmp = (PTM *)malloc(sizeof(PTM));
    if (!ptm_tmp) {
        printf("Error allocating memory for ptm array!\n");
        exit(1);
    }

    // Read modifications into ptm array from cnfg file.
    while (fgets(modStr, sizeof(modStr), cnfg) != NULL) {

        // Change code to read from new type config file
        if (strstr(modStr, "_varMod")) {
            PTM *tmp = realloc(ptm_tmp, (1 + i)*sizeof(PTM));
            if (!tmp) {
                printf("Error reallocating memory for ptm array!\n");
                free(ptm_tmp);
                exit(1);
            }
            ptm_tmp = tmp;

            sscanf(modStr, "_varMod\t%s\t%lf\t%s", ptm_tmp[i].name, &ptm_tmp[i].delta_mass, ptm_tmp[i].site);
            i++;
        } 
        
    }

    double model_params[4];
    // if not first search, recalibrate masses. This suposes that first search has been run and there is corresponding file
    if (first_search == 2) {
        for (int j = 0; j < 4; j++) {
            model_params[j] = 0;
        }
        char str44[128] = {};
        strcat(str44, argv[1]);
        strcat(str44, ".1st");
        strcat(str44, ".txt");
        mass_recalibration(i, str44, model_params, min_score, max_data);
    }
    
    // for debugging
    printf("precTol = %lf\n", precTol);
    printf("tol = %lf\n", tol);
    printf("maxMod = %d\n", maxMod);
    printf("FASTA files are:\n");
    for (int j = 0; j < fasta_cntr; j++) printf("%s\n", fasta[j]);
    printf("Variable modifications are:\n");
    for (int j = 0; j < i; j++) printf("%s\t%lf\t%s\n", ptm_tmp[j].name, ptm_tmp[j].delta_mass, ptm_tmp[j].site);
    //return 0;

    // Create counts[i] array to hold number of modifications of each type..
    int counts[i];
    for (int w = 0; w < i; w++) counts[w] = 0;

    // Open sqlite database.
    sqlite3 *db;
    sqlite3_stmt *ppStmt;
    int rc;

    //printf("Finished copying to in-memory database\n");

    FILE *mms_ptr;
    double lower;
    double upper;
    //int ct = 0;
    int ct1 = 0;
    printf("Processing %s...\n", argv[1]);
    
    // Open results file for writing.
    /*char resFile[128] = {};
    strcat(resFile, argv[1]);
    if (first_search) strcat(resFile, ".1st");
    strcat(resFile, ".txt");*/
    rc = sqlite3_open(argv[2], &db);
    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        exit(1);
    }
    
    // copy to in memory database
    // now need to copy some data to an in-memory db that will be used in eval_obj()
    sqlite3 *db1;
    sqlite3_backup *pBackup;

    rc = sqlite3_open(":memory:", &db1); // need to write error handling code later
    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db1));
        sqlite3_close(db1);
        exit(1);
    }
    
    pBackup = sqlite3_backup_init(db1, "main", db, "main");
    if (pBackup) {
	    (void)sqlite3_backup_step(pBackup, -1);
	    (void)sqlite3_backup_finish(pBackup);
	    printf("Finished copying to in-memory database\n");
    }
    else {
        rc = sqlite3_errcode(db1);
        printf("Error rc = %d\n", rc);
	exit(1);
    }
    //if (sqlite3_errmsg(db) sqlite3_free((void *)sqlite3_errmsg(db));
    sqlite3_close(db);

    // debug
    /*sqlite3_close(db1);
    free(ptm_tmp);
    xit(1);*/

    // Read data struct by struct, search and write results to file.
    // 09.12.2023: now refactoring to write results at the end
    if ((mms_ptr = fopen(argv[1], "rb")) == NULL) {
        printf("Error opening binary file for reading.\n");
        exit(1);
    }
    // update 06.09.2021: will count spectra to determine N_spectra here
    spectrum spec1;
    while (fread(&spec1, sizeof(spectrum), 1, mms_ptr) == 1) {
        N_spectra++;
    }

    // 9.01.2024: starting the first search implementation
    /*
     * Need to:
     * 1. Read first search config parameters—change the parsing function to do this
     * 2. Run a first search but do not write results to file, collect them in arrays
     *    could be fixed length arrays
     * 3. least_squares to fit the model
     * 4. Change the scorer to correct PrecMasses according to the model and run as before
     *
     * Will have to refactor eventually. Should make a scorer function to be called with 
     * parameter struct so I can reuse it for first and main search.
     * So, 
     * 1. a scorer function to generate results array. Should take parameter to tell it if it is
     *    first or main search.
     * 2. a least_square function to fit model
     * 3. New main will just call scorer and print results to file
     */
    
     

    // 9.12.2023: now declare results array but allocate it on the heap beacause it will be too large for the stack
    // don't forget to free ti!!
    results *results = NULL;
    //results results[2*N_spectra];

    if ((res = fopen(resFile, "w")) == NULL) {
        printf("Error opening results file for writing.\n");
        exit(1);
    }

    fprintf(res, "Total\tMatched\tIntensity_matched\tq"); // Need matched and total ion numbers (18.06.2020).
    for (int ctr = 0; ctr < i; ctr++) fprintf(res, "\t%s probs", ptm_tmp[ctr].name); // later should put results in data struct perhaps
    fprintf(res, "\tFragments\tScan\tPrecScan\tRetTime\tSequence\tCharge\tPrecMass\tMassError\tScore\tProteins");
    for (int ctr = 0; ctr < i; ctr++) fprintf(res, "\t%s", ptm_tmp[ctr].name);
    fprintf(res, "\tnMod\tMissed\tDecoy\tContam\tn_seq\tscore_threshold\n");

    //spectrum spec1;
    char sql_inst[1024]; 
    int specNo;

    // starting omp parallelization
    // update 06.09.2021: will go back to for loop as tasks are hard to get to work
    
    #pragma omp parallel default(none) shared(score_table, model_params, ptm_tmp, mms_ptr, maxMod, precTol, tol, precTol_1st, first_search, check_alt_iso, db1, i, res, N_spectra, ct1) private(specNo, spec1, sql_inst, lower, upper, counts, seq, ppStmt, rc, proteins, str3, item, item1, mod, seq_count) firstprivate(results1, results2, results, res_count, str22)

    { 
    
    // need to malloc ptm, item and item1
    // 9.12.2023: also arrays to hold results
    #pragma omp critical
    {
    // Create item pointers and allocate memory.
        item = malloc(sizeof(*item));  // This might need to go to scorer_mod.c and be passed as argument.
        if (!item) {
            printf("Error allocating memory for item!\n");
            exit(1);
        }

        // for debugging
        //printf("Created item in thread %d address is %p\n", omp_get_thread_num(), (void *)item);
 
        item1 = malloc(sizeof(*item1)); // This might need to go to scorer_mod.c and be passed as argument.
        if (!item1) {
            printf("Error allocating memory for item!\n");
            exit(1);
        }

        ptm = (PTM *) malloc(i*sizeof(PTM));
        if (!ptm) {
            printf("Error allocating memory for ptm array!\n");
            exit(1);
        }
        for (int j = 0; j < i; j++) {
            strcpy(ptm[j].name, ptm_tmp[j].name);
            strcpy(ptm[j].site, ptm_tmp[j].site);
            ptm[j].delta_mass = ptm_tmp[j].delta_mass;
        }
        // for debugging
        //printf("Created ptm in thread %d address is %p\n", omp_get_thread_num(), (void *)ptm);
        //printf("Name of ptm[0] is %s\n", ptm[0].name);
        //printf("delta mass of ptm[0] is %lf\n", ptm[0].delta_mass);

        /*mod = (mod_ind *)malloc(sizeof(*mod));
        if (!mod) {
            printf("Error allocating memory for mod array!\n");
            exit(1);
        }*/
        mod = NULL;
        //printf("Created mod array in thread %d address is %p\n", omp_get_thread_num(), (void *)mod);
        results = malloc(2*N_spectra*sizeof(*results));
    }


    #pragma omp for reduction(+:ct1)
        for (specNo = 0; specNo < N_spectra; specNo++) {
            #pragma omp critical
            {
                fseek(mms_ptr, sizeof(spectrum) * specNo, SEEK_SET);
                fread(&spec1, sizeof(spectrum), 1, mms_ptr);
            }
            //lower = spec1.precMass - (precTol*spec1.precMass)/1000000;
            //upper = spec1.precMass + (precTol*spec1.precMass)/1000000;

	    /*
	     * Update 02.12.2023: will implement alternative precursor isotopes since low-intensity spectra seems to not be 
	     * handled well. MaxQuant assigns more spectra than scorer_omp as of present. But scorer_omp assigns more or same
	     * scribble spectra in pull down data. This suggests that abundant peptide ions are handled well, presumably because
	     * isotope patterns are good and precursor mass is easy to determine. Obviously, it will slow down execution but is 
	     * worth trying.
	     */
            if (first_search == 1) {
		lower = spec1.precMass - (precTol_1st*spec1.precMass)/1000000;
                upper = spec1.precMass + (precTol_1st*spec1.precMass)/1000000;
	        sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf) and missed = 0 and mod_pep.mass = peptides.mass", lower, upper);
	    }
	    else if (first_search == 2) {
                // redefine lower and upper using the model
                double new_mass = spec1.precMass - spec1.precMass*(model_params[0] + model_params[1]*spec1.precMass + 
				model_params[2]*spec1.retTime + model_params[3]*spec1.precMass*spec1.retTime)/1000000;
                lower = new_mass - (precTol*new_mass)/1000000;
		upper = new_mass + (precTol*new_mass)/1000000;
		if (check_alt_iso == 0) {
                    sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf)", lower, upper);
                }
                else {
                    sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf) or (mod_pep.mass between %lf and %lf)", lower, upper, lower + Pr, upper + Pr);
                }
	    }
            else {	    
                lower = spec1.precMass - (precTol*spec1.precMass)/1000000;
                upper = spec1.precMass + (precTol*spec1.precMass)/1000000;
		if (check_alt_iso == 0) {
		    sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf)", lower, upper);
		}
		else {
                    sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf) or (mod_pep.mass between %lf and %lf)", lower, upper, lower + Pr, upper + Pr);
		}
            }
	    /*if (check_alt_iso == 0) {
                sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf)", lower, upper);
            }
	    else if (check_alt_iso == 2) sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf) or (mod_pep.mass between %lf and %lf) or (mod_pep.mass between %lf and %lf)  or (mod_pep.mass between %lf and %lf)  or (mod_pep.mass between %lf and %lf)", lower, upper, lower - Pr, upper - Pr, lower - Pr - Pr, upper - Pr - Pr, lower + Pr, upper + Pr, lower + Pr + Pr, upper + Pr + Pr);
	    else sprintf(sql_inst, "select * from mod_pep left join peptides on mod_pep.pep_id=peptides.rowid where (mod_pep.mass between %lf and %lf) or (mod_pep.mass between %lf and %lf) or (mod_pep.mass between %lf and %lf)", lower, upper, lower - Pr, upper - Pr, lower + Pr, upper + Pr);*/

            //printf("Searching for %lf...\n", spec1.precMass); 
            rc = sqlite3_prepare_v2(db1, sql_inst, -1, &ppStmt, 0);

            while ((rc = sqlite3_step(ppStmt)) != 101) {
                if (rc == 100) {
                    int rowid = sqlite3_column_int(ppStmt, 0);
                    for (int ctr = 0; ctr < i; ctr++) counts[ctr] = sqlite3_column_int(ppStmt, ctr + 2);
                    for (int ctr = 0; ctr < i; ctr++) ptm[ctr].count = sqlite3_column_int(ppStmt, ctr + 2);
                    strcpy(seq, (char *)sqlite3_column_text(ppStmt, 2 + i));
                    strcpy(proteins, (char *)sqlite3_column_text(ppStmt, 4 + i));
                    int missed = sqlite3_column_int(ppStmt, 5 + i);
                    int decoy = sqlite3_column_int(ppStmt, 6 + i);
                    int contam = sqlite3_column_int(ppStmt, 7 + i);
                    int count_sum = 0; // 9.12.2023: this seems not necessary. mod_pep is lready generated to comply with maxMod
                    for (int ctr = 0; ctr < i; ctr++) {
                        count_sum += counts[ctr];
                        //printf("%s is %d\n", ptm[ctr].name, ptm[ctr].count);
                    }

                    if (count_sum <= maxMod) { 
                        bst(counts, seq, 42, item, item1, mod, i);
			callback(rowid, counts, i, seq, &spec1, proteins, item, item1, mod, missed, 
					decoy, contam, count_sum, &results1, &results2, str3, score_table);
			STdestroy(head);
			free(z); //16.12.2023: took me few days to arrive at this, now valrgind is happy!
                    }
		    seq_count++;
                }
                    
                else {
                    printf("Error executing sql statement: %d.\n", rc);
                    exit(1);
                }
                

            }
            ct1++;
            sqlite3_finalize(ppStmt);
           
	    // 9.12.2023:now write results to results array
	    //printf("res_count = %d\n", res_count);

            sprintf(results[res_count].data, "%s", results1.data); 
	    if (strcmp(results1.data, "")) {
                sprintf(str22, "%d\t", seq_count);
		strcat(results[res_count].data, str22);
		sprintf(str22, "%lf\n", -10*log10(0.01/seq_count));
		strcat(results[res_count].data, str22);
	    }
	    res_count++;
            sprintf(results[res_count].data, "%s", results2.data);
            if (strcmp(results2.data, "")) { // there wa a bug here. instead of results2.data, results1.data was copied from above
                sprintf(str22, "%d\t", seq_count);
		strcat(results[res_count].data, str22);
                sprintf(str22, "%lf\n", -10*log10(0.01/seq_count));
		strcat(results[res_count].data, str22);
         
	    }

	    res_count++;

            // set results.score to 0 and wipe clean data fiels
            results1.score = 0;
            results2.score = 0;
            strcpy(results1.data, ""); // necessary otherwise scans that do not have hits are replaced by previous scans data
            strcpy(results2.data, "");
	    seq_count = 0;

        }
	#pragma omp critical
        {
            // 10.12.2023: now print results to file
            //printf("Finished searching in this thread! Now printing results to file...\n");
            for (int i = 0; i < res_count; i++) {
                fprintf(res, "%s", results[i].data);
            }
        }        
	    free(results);
            free(ptm);
            free(item);
            free(item1);
            free(mod);
            mpfr_free_cache();
    }


    (void)sqlite3_close(db1);
    fclose(mms_ptr);
    fclose(cnfg);
    fclose(res);
    free(ptm_tmp);  
    printf("Processed %d MS/MS spectra.\n", ct1);
    
    return 0;
}

