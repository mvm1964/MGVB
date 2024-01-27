/* Bayesian Network for computing modification site probabilities

   Metodi V. Metodiev
   University of Essex, United Kingdom


   Update 03.08.2022: the ideas described above are interesting but combining fragments is
   not correct with respect to robability theory as pairs are not independent evidence.
   Will change now to the model testing approach as developed in July 2022. Will use again a BN 
   but much sompler. Each fragment will be assessed whether it is compatible with a secific 
   configuration (hypothesis, model). See python code for more on the algorithm.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "pepMass.h"
#include "item_open.h"
#include "compute_site_prob_optimized.h"

/*typedef struct Results {
    double score;
    char data[4096];
} results;*/

int get_n(char *seq, int i) {
    int str_size = strlen(seq);
    int n = 0;
    for (int j = 0; j < str_size; j++) {
        if (strchr(ptm[i].site, seq[j]) != NULL) n++;
    } 
    return n;
}

int get_m(char *enc, int i) {
    int str_size = strlen(enc);
    int m = 0;
    for (int j = 0; j < str_size; j++) {
        if (enc[j] == i + 48) m++;
    }
    return m;
}

// generate model: recursive function to generate model specifications
void model_spec(int m, int n, int j, int *model, int specs[][n], int k, int *ctr) {
    /*
       m: number of modifications of particular type on peptide
       n: number of candidate sites
       j: start of sequence for recursion
       model: binary array with sites ids
       specs: 2d array with models
       k: counter for modifications
       ctr: counter for models
    */

    if (k == m) {
        // set all entries to 0
        for (int w = 0; w < n; w++) specs[*ctr][w] = 0;

        // assign 1 to specs[][model[w]]. Model has m entries holding the index values (i)
        for (int w = 0; w < m; w++) specs[*ctr][model[w]] = 1;

        //for (int w = 0; w < n; w++) printf("%d ", specs[*ctr][w]);
        //printf("\n");
        //printf("ctr = %d\n", *ctr);
        (*ctr)++;
        return;
    }

    for (int i = j; i < n; i++) {
        model[k++] = i;
        model_spec(m, n, i + 1, model, specs, k, ctr);
        k--;

    }
}

// test_compat function: test if fragment from item array is compatible with a model
int test_compat(int *model, Item item, int nmod, int nsite, int mod_type) {
    /*
       model: binary array of modification sites
       item: fragment with information about type, idex, modifications. Modifications
             are in the enc field of the struct
       nmod: number of modifications on peptide
       nsite: numebr of candidate sites on peptide
    */
    
    // need to determine nmod, the number of modifications 

    int test = 0;
    if (strchr(item->fragType, 'b')) {
        int predicted_m = 0;        
        for (int i = 0; i < nsite; i++) predicted_m += model[i];
        if (predicted_m == nmod) test = 1;
    }
    else if (strchr(item->fragType, 'y')) {
        int predicted_m = 0;
        for (int i = nsite - get_n(item->seq, mod_type); i < nsite; i++) predicted_m += model[i];
        if (predicted_m == nmod) test = 1;
    }
    return test;
}

// main function
void compute_site_prob_optimized(results *results_any, char *str3, int mod_type, int counts, char *seq, Item *item_array, int size) {
    
    int m = counts;
    int n = 0;
    // Now determine n using a function and the seq
    n = get_n(seq, mod_type);
    // need to test if m == n and if so return with unit vector of length m for probs
    if (m == n) {
        for (int i = 0; i < n; i++) {
            sprintf(str3, "%lf,", 1.0);
            strcat(results_any->data, str3);
        }
        strcat(results_any->data, "\t");
        return;        
    }
    // update 10.08.2022: perhaps should check if n > some number and only then allocate on the heap.
    // this should make it faster    
    // now generate the models. Will follow the python notebook but later will use bitmasks for the combinations
    // first need to compute combinations number
    mpz_t comb;
    int comb1;
    mpz_init_set_ui(comb, 0);
    mpz_bin_uiui(comb, n, m);
    //mpf_set_z(comb1, comb);
    comb1 = mpz_get_d(comb);

    // for debugging
    //printf("comb1 = %d\n", comb1);

    mpz_clear(comb);
    
    // need to allocate this array on the heap because stack cannot accomodate more than 34220 combs
    //int models[comb1][n];
    //int *models = NULL;
    /*if (n > 20) {
        (*models)[n] = calloc(comb1, sizeof *models);
    }
    else {
        models[comb1][n];
    }*/
    int (*models)[n] = calloc(comb1, sizeof *models);
    int frag_m = 0;
    
    // variable to hold the likelihood
    mpfr_t L[comb1];
    double tmp = (double)(1.0/comb1); // very important to use 1.0 not 1
    for (int i = 0; i < comb1; i++) {
        mpfr_init(L[i]);
        //printf("tmp = %lf\n", tmp);
        mpfr_set_d(L[i], tmp, MPFR_RNDN);
    }

    // for debugging
    //for (int i = 0; i < comb1; i++) printf("L is %lf\n", mpfr_get_d(L[i], MPFR_RNDN));

    int model[n];
    for (int i = 0; i < n; i++) model[i] = 0; 
    int ctr = 0;
    // generate model matrix
    model_spec( m, n, 0, model, models, 0, &ctr);    
    
    // now loop over item_array and update site probabilities 
    for (int i = 0; i < size; i++){
        frag_m = get_m(item_array[i]->enc, mod_type);
        for (int j = 0; j < comb1; j++) {
            if (test_compat(models[j], item_array[i], frag_m, n, mod_type)) {
                mpfr_mul_d(L[j], L[j], 0.8, MPFR_RNDN); // probability should be settable in config files
            }
            else mpfr_mul_d(L[j], L[j], 0.15, MPFR_RNDN);
        }    
    }

    mpfr_t sum_w;
    mpfr_init_set_d(sum_w, 0, MPFR_RNDN);
    for (int i = 0; i < comb1; i++) mpfr_add(sum_w, sum_w, L[i],  MPFR_RNDN);

    // now normalise and write results to file and return
    for (int w = 0; w < comb1; w++) {
        mpfr_div(L[w], L[w], sum_w, MPFR_RNDN);              

    }

    // compute site probabilities
    double prob_tmp;
    double prob[n];
    for (int i = 0; i < n; i++) prob[i] = 0;
    for (int i = 0; i < comb1; i++) {
        prob_tmp = mpfr_get_d(L[i], MPFR_RNDN);
        for (int j = 0; j < n; j++) {
            if (models[i][j] == 1) prob[j] += prob_tmp;
        }
    }

    for (int i = 0; i < n; i++) {
        sprintf(str3, "%lf,", prob[i]);
        strcat(results_any->data, str3);
    }

    // for debugging
    //printf("%s\n", str3);
 
    strcat(results_any->data, "\t");    

    mpfr_clear(sum_w);
    for (int i = 0; i < comb1; i++) mpfr_clear(L[i]);
    //if (n > 20) free(models);
    free(models);

    return;
   
}
