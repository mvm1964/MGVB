/*
  BST functions
  
  Metodi V. Metodiev, 2020, UK

  Generates combinations of possible modified peptides.
  Uses the data returned from db to generate combinations:
  
    - if there are 2 phospho, one M_Ox and 1 N-term acetylation
      this is coded  as {2, 1, 1}. 
    - This is then used together with the peptide sequence to 
      generate the combinations.
    - Uses PTM and mod_ind types defined in pepMass.h.

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pepMass.h" // changed on 11.09.2022
#include "item_open.h"
#include "bst2023.h"
#include <math.h>
#include <omp.h>

link head, z;
int comb_index = 0;
#pragma omp threadprivate(head, z, comb_index)
 
// function that grows flexible int arrays inside items 
 Item grow_array(Item array) {
    int new_max = array->max_flex + 100;
    Item tmp = realloc(array, sizeof(*array) + new_max*sizeof(int));
    if (!tmp) { // Should consider calling a custom exit function to avoid memory leaks on exit (22.06.2020).
        printf("Error reallocating memory for comb_ind array!\n");
	printf("comb_index, seq, enc are %d, %s, %s\n", comb_index, array->seq, array->enc);
	free(array->key);
	//free(array->comb_ind);		
        free(array);
        exit(1);
	//return(NULL);
    }
    array = tmp;
    array->max_flex = new_max;
    for (int i = array->flex_count; i < array->max_flex; i++) array->comb_ind[i] = 0;

    // for debuging
    //printf("in grow_array() comb_index, max_flex, and flex_count are: %d, %d, %d\n", comb_index, array->max_flex, array->flex_count);
    //for (int i = 0; i< array->flex_count; i++) printf("comb_ind[%d] = %d\n", i, array->comb_ind[i]);
    return array;
}	


/* Recursive function to generate C(k, n) of an modifications array. Uses a graph traversal algorithm.
   Update 20.06.2020: Need to handle unmodified peptides. Now it crashes because k is 0.

   Update 17.08.2021:
   Add head and z links to the argument lists of bst, comb_kN, mid_xcorr, STinit,
     STcount, STsearch, STsearch_ins,  STsearch, int_STsearch, STinsert
*/
void comb_kN(mod_ind *arr, int *ind_array, int *data, int k, int n, char *str1, 
                    int start, int ind, int i, char *seq, int seq_len, Item item, Item item1) {

    /* An array to store the modifications counts. It has the same number of elements 
       as there are modifications. NB: Consider defining it as external to speed up execution.
    */

    int counts[i];
    for (int j = 0; j < i; j++) counts[j] = 0;

    // For debugging
    //for (int p = 0; p < i; p++) printf("%s count is %d\n", ptm[p].name, ptm[p].count);
    //exit(1);

    // Indicator var to tell if string is productive.
    int indicator = 0;

    // If combination is ready it is processed.
    if (ind == k) {
        
        // Test if modification set meets requirements.
        for (int j = 0; j < k; j++) {
            for (int t = 0; t < i; t++){ 
                if (str1[data[j]] == t + 48) { //&& counts[t] < ptm[t].count) { // 12.07.2020: not sure if this is correct. Should probably remove the && condition. 
                    counts[t]++;
                }
            }
        }


        for (int j = 0; j < i; j++) {
            if (counts[j] != ptm[j].count) indicator++;
        }
        if (indicator == 0) {
            mid_xcor(arr, data, str1, seq, k, item, item1);
        }
        
        return;
    }

    // From http://exceptional-code.blogspot.com/2012/09/generating-all-permutations.html
    for (int j = start; j < n; ++j) {
        data[ind++] = ind_array[j];
        // Recur with str + 1.
        comb_kN(arr, ind_array, data, k, n, str1, j + 1, ind, i, seq, seq_len, item, item1);
        --ind; // Is this necessary. Can't we just recur with ++ind.
    }
}

// This is to add to BST and print modified strings.
void mid_xcor(mod_ind arr[], int *data, char *str1, char *seq, int r, Item item, Item item1) {
    int str_len = strlen(seq);
    char str2[1024];
    strcpy(str2, str1);
    for (int j = 0; j < str_len; j++) str2[j] = '.';
    for (int v = 0; v < r; v++) str2[data[v]] = str1[data[v]];

    // Populate the item fields and build BST.
    double extraMass = 0;
    strcpy(item->seq, seq);
    strcpy(item->enc, str2);
    for (int j = 0; j < str_len; j++) {
        if (str2[j] != '.') extraMass = extraMass + ptm[str2[j] - 48].delta_mass;
    }

    double Mass = pepMass(seq) + extraMass;
    item->key = &Mass; 

    // for debuging
    //printf("encoding of this combination is %s\nsequence is %s\n", item->enc, item->seq);

    bst_y(head, item, 0, str_len- 1, item1);  // it has to be strlen(seq) - 1
    // for debugging
    //printf("Added y fragments. Sequence is %s\n", seq);

    // for debuging will not do b fragmets for now
    bst_b(head, item, 0, str_len - 1, item1);
    // for debugging
    //printf("Added b fragments\n");
    //printf("n4 = %d\n", head->N);    

    // 30.10.2022: if I want individual searches, need to call a scoring function here
    // Some kind of lightwait xcorr that will use the tree to calculate score and then destroy it.

    // increment comb_index
    comb_index++;
    
    // for debuging
    //printf("encoding of this combination is %s\nsequence is %s\n", item->enc, item->seq);

}

/*
  Main function: it builds the bst and returns a pointer to it (no, it just builds it).
  It receives as arguments a counts array, number of elements in 
  the array, a peptide sequence string, and a PTM array.
  Update 06.06.2020: consider moving this within scorer_mod.c, could
  simplify implementation for now. Later could make it into a function. 
*/
void bst(int *counts, char *seq, int n, Item item, Item item1, mod_ind *mod, int i) {

    // 08.12.2023: n does not seem to be used at all???
    
    // for debugging
    // printf("in bst() comb_index is %d\n", comb_index);
    //comb_index = 1;
    //int i = 0; this comes as function argument now (07.06.2020).
    //char *token = NULL; 
    //char str[1024]; // Container for reading modfication string from file.
    char seq1[1024]; // Container for peptide sequence.

    int start = 0;
    int ind = 0;

    // Create item pointers and allocate memory.
    /*Item item = malloc(sizeof(*item));  // This might need to go to scorer_mod.c and be passed as argument.
    if (!item) {
        printf("Error allocating memory for item!\n");
        exit(1);
    }

    Item item1 = malloc(sizeof(*item1)); // This might need to go to scorer_mod.c and be passed as argument.
    if (!item1) {
        printf("Error allocating memory for item!\n");
        exit(1);
    }*/
     
    // Create ptm array and populate it from config file, assuming file is modifications.txt.
    // Update 14.05.2020: this should be moved to scorer_mod.c.
    /*PTM *ptm = malloc(sizeof(*ptm));
    if (!ptm) {
        printf("Error allocating memory for ptm array!\n");
        exit(1);
    }*/

    // Array to hold modifications indices. Will realloc when data comes from sqlite3 db.
    // Update 14.05.2020: this should be moved to scorer_mod.c.
    /*mod_ind *mod = malloc(sizeof(*mod));
    if (!mod) {
        printf("Error allocating memory for mod array!\n");
        exit(1);
    }*/

    /*FILE *cnfg; // Update 14.05.2020: this should be moved to scorer_mod.c.
    if ((cnfg = fopen("modifications.txt", "r")) == NULL) {
        printf("Error openning modifications file!\n");
        exit(1);
    }*/

    // Update 14.05.2020: this should be moved to scorer_mod.c.
    /*while (fgets(str, sizeof(str), cnfg) != NULL) {
        PTM *tmp = realloc(ptm,  (i + 1)*sizeof(PTM));
        if (!tmp) {
            printf("Error reallocating memory for ptm array!\n");
            free(ptm);
            exit(1);
        }
        ptm = tmp;
        token = strtok(str, "\t");
        strcpy(ptm[i].name, token);
        token = strtok(NULL, "\t");
        ptm[i].delta_mass = atof(token);
        token = strtok(NULL, "\t");
        strcpy(ptm[i].site, token);
        ptm[i].count = 0;
        //strtok(NULL, "\t");
        i++;
    }
    printf("There are %d modifications to search for.\n", i); // Update 14.05.2020: this should be moved to scorer_mod.c.
    */
   
    // Read sequence and find modification indices. 

    //int counts[] = {2,1,1}; // In real app this will come from db.
    

    // Update counts field of ptm structs. // Update 14.05.2020: this should be moved to scorer_mod.c.
    /*for (int j = 0; j < i; j++) {
        ptm[j].count = counts[j];
    }*/

    //strcpy(seq, argv[2]); // Update 14.05.2020: seq should be a function argument now..
    int seq_len = strlen(seq);
    // printf("Length is %d\n", seq_len);
    strcpy(seq1, seq);  

    //int ctr = 0;
    for (int j = 0; j < seq_len; j++) seq1[j] = '.';
    //printf("%s\n", seq1);

    /*
       Update 14.08.2021: need to change it to account for overlapping specificities although this will slow down 
       and is only worth doing for the open search. Should generate combinations of encoding strings and use all of
       them to generate the tree. Perhaps it would be only needed for now to change the mod_comb data to contain
       only combinations of non-overlapping modifications.
    */
   
    // Generate the encoding string. 
    for (int v = 0; v < seq_len; v++) {
        for (int j = 0; j < i; j++) {
            if (strchr(ptm[j].site, seq[v]) != NULL) {
                seq1[v] = j + 48; 
                
                            
            }
        }
    }
    
    // For debugging
    //for (int ctr = 0; ctr < i; ctr++) printf("%s is %d\n", ptm[ctr].name, counts[ctr]);

    // Now call comb_kN.
    //printf("Sending %s with lenght %lu\n", seq, strlen(seq));
    
    // Find total number of modified residues.
    // Need to modify to handle unmodified peptides. Should I use the original xcorr function?
    int k = 0; // Now this is determined in scorer_mod. Should modify code to receive it from there to speed up.
    for (int j = 0; j < i; j++) k = k + counts[j];     
    //printf("%d residues are modified.\n", k);
    
    // For debugging only
    //if (k == 0) return; 

    // Find total number of modifiable residues (20.06.2020: if k = 0 then should go directly to xcorr but should place 
    // this test in scorer_mod before calling specific scoring functions).
    int j1 = 0;
    for (int j = 0; j < seq_len; j++) {
        if (seq1[j] != '.') {
            j1++;
        }

    }
    

    // Reallocate mod array memory.
    int j = 0; 
    while (j < j1) {
        mod_ind *tmp_mod = realloc(mod,  (j + 1)*sizeof(mod_ind));
        if (!tmp_mod) { // Should consider calling a custom exit function to avoid memory leaks on exit (22.06.2020).
            printf("Error reallocating memory for mod array!\n");
            free(mod);
            exit(1);
        }
        mod = tmp_mod; 
        //printf("%d\n", j);
        j++;
    }
    // Update mod array.
    int j2 = 0;
    for (int j = 0; j < seq_len; j++) {
        if (seq1[j] != '.') {
            mod[j2].index = j;
            mod[j2].type = seq1[j];
            j2++;
        }
    }
    
    // Build the index array.
    int ind_array[j1];
    for (int j = 0; j < j1; j++) ind_array[j] = mod[j].index;

    //for (int j = 0; j < j1; j++) printf("%d\n", ind_array[j]);
    //printf("\n");

    int data[k]; // Array to store combinations.
    
    //for (int p = 0; p < k; p++) data[p] = 0; // for debugging

    // For debugging
    //for (int p = 0; p < i; p++) printf("%s count is %d\n", ptm[p].name, ptm[p].count);

    // set comb_index to 0
    comb_index = 0;
    
    // Create BST
    STinit();
    //printf("BST created.\n"); // For debugging
    //int start = 0;
    //int ind = 0; 
    comb_kN(mod, ind_array, data, k, j1, seq1, start, ind, i, seq, seq_len, item, item1);

    /* 
    // Test BST.
    for (int j = 1; j < 5; j++) {
        double mass = j;
        printf("Inserting %s\n", argv[2] + j);
        STinsert(newItem(j, argv[2] + j, "...", &mass));
    }
    */

    //printf("The binary search tree has %d nodes.\n", head->N);
    //printf("Head sequence is %s\n", head->item->seq);
    //STprint(head, z);
    /*
    //Some tests of the interval search.
    char resp;
    double x1 = 1015.290;
    double x2 = 1015.490540;
    printf("Would you care for some tests?\n");
    printf("Enter low: ");
    scanf("%lf", &x1);
    printf("Enter high: ");
    scanf("%lf", &x2);
    //Item result = int_STsearch(x1, x2);
    rangeSelect(head, x1, x2);
    */
    /*
    if (result != NULLitem) {
        printf("%s\n%s\n%lf\n", result->seq, result->enc, *(result->key));
    }
    
    else printf("Not found in BST.\n");
    */

    /*fclose(cnfg); // This should go to scorer_mod.c
    free(mod);    // As above
    free(ptm);    // Not sure if it should be freed and recreated after each sequence from db.  
    free(item);   // As above
    free(item1);  // As above
    STdestroy(head); // Should go to scoring function.
    free(z);*/  
    free(mod); // It has to be freed here to prevent leaks. For whatever reason...
    //free(z);  // Same as mod. Update 15.08.2022: now it is in STDestroy()
    // for debugging
    //printf("in bst() comb_index is %d\n", comb_index);
    
    return;

}

// Functions definitions.

// Creates a new tree node.
link NEW(Item item, link l, link r, int N) {
    link x = malloc(sizeof *x);
    x->item = item;
    x->l = l;
    x->r = r;
    x->N = N;
    return x;
}

// Initializes the tree.
void STinit() {
    head = (z = NEW(NULLitem, NULL, NULL, 0));
}

// Counts nodes in subtree.
int STcount(link head) {
    return head->N;
}

// Recursive search for mass.
Item searchR(link h,  Key v) {
    //Key t = key(h->item);
    if (h == z) return NULLitem;
    Key t = key(h->item);
    if (eq(v, t)) return h->item;
    if (less(v, t)) return searchR(h->l, v);
    else return searchR(h->r, v); 
}

// Modified searchR for inserting: considers enc string as well as key and seq.
Item searchR_ins(link h, Key v, char *seq, char *enc) { // seq is actually enc
    //Key t = key(h->item);
    if (h == z) return NULLitem;
    Key t = key(h->item);
    /* Update 18.07.2020: need to cover cases where different modifications sets have the same delta mass
       should write a new Item function that compares keys, seq and modifications to decide if the new item 
       should be added. Can use the code from bst_y(), which scans the enc string and collect modification info
    */

    // For debugging
    //printf("seq = %s\n", seq);
    //printf("enc = %s\n", enc);
    //printf("h->item->seq = %s\n", h->item->seq);


    if (eq(v, t) && !strcmp(h->item->seq, seq) && !item_comp(h->item, enc)) {
        // need to reallocate h->item->comb_ind
	if (h->item->flex_count >= h->item->max_flex) {
	    Item tmp = NULL;
            tmp = grow_array(h->item);
	    h->item = tmp;
	}
        return h->item;
    }

    /*if (eq(v, t) && item_comp(h->item, seq)) { // This part slows down the execution a lot
        searchR_ins(h->l, v, seq);               // because it is wrong. It inserts the frag no mater if it exists on the
        searchR_ins(h->r, v, seq);               // other branch. That is why it uses double the memory. 
    }*/
    
    /* For debugging: it does not find any with maxMod 6 and phos, acetyl and metOx. Should check from time to time
       with other modifications and if never finds any should probably drop the item_com call in searchR_ins().
    if (eq(v, t) && !strcmp(h->item->seq, seq) && item_comp(h->item, enc)) {
        printf("Found it!\n");
        printf("New seq is\t%s\n", seq);
        printf("New enc is\t%s\n", enc); 
        printf("h->item->seq is\t%s\n", h->item->seq);
        printf("h->item->enc is\t%s\n", h->item->enc);
    }
    */

    if (less(v, t)) return searchR_ins(h->l, v, seq, enc);           
    else return searchR_ins(h->r, v, seq, enc);
}

// STsearch midware.
Item STsearch(Key v) {
    return searchR(head, v);
}

// Modified STsearch for inserting: considers seq string as well as key.
Item STsearch_ins(Key v, char *seq, char *enc) {
    return searchR_ins(head, v, seq, enc);
}

// Interval search functions. Update 12.07.2020: this does not seem right. Need to traverse and add all matching leaves.
Item int_searchR(link h, double low, double high) {
    if (h == z) return NULLitem;
    Key t = key(h->item);
    if (int_eq(low, high, t)) return h->item; // Update 12.07.2020: instead of returning it should add the item to an array and recur.
    if (int_less(t, low)) return int_searchR(h->r, low, high);
    else return int_searchR(h->l, low, high);
}

/* Update 12.07.2020: New interval search function. It should save all matching fragments to a new array of items 
   For each BST there should be one item_array containing all fragments that match the spectrum masses. 
   The function receives the head link, the low and high mass values, the pointer to the item_array and the current size
   of the item_array. It then searches for matching items in the BST and adds them to the item_array after reallocating it.
   Question: what happens with n_item in the recursive calls?

   Update 16.07.2020: The in order traversal is not efficient for interval_search() because it visits all leaves.
   Need to change the algorithm to ignore branches if they are outside the interval:
   1. Check for base condition and return if met.
   2. Check if root is in range. If so, recur to left and right.
   3. If root less than low recur to right.
   4. If root greater than high recur to left.

   Update 11.09.2022: implementing H2O and NH3 loss, and 2+ fragments search
   Need to pass frmz pointer and the index of the matched fragment mass so they
   can be used in the functions that search for H2O and NH3 loss
*/
int interval_search(link h, Item **item_array, int *n_item, double *frmz, 
		int dblchrg, double tol, int *y_fr, int *b_fr, double *frint, int *n2, int *n4, 
		int *n2_array, int numMod, int *n3_array, int *n4_array, int *n5_array, int n_item1, Item *item_array1, int n) {
    
    /*
     * Start the new implementation 01.01.2024
     *
     * 1. begin with item_array1[0] and frmz[0]: set index to 0
     * 2. if match, process the item as in old implementation and update the index variable 
     *    if no match proceed to next element of item_array1 and so on
     *
     */    

    int ctr22 = 0;
    int ctr33 = 0;
    int start_index = 0;
    int match = 0;
    char str4[1024];

    // debug
    //printf("frmz is %d elements long\n", n);
    
    if (dblchrg) {
        for (int i = 0; i < n; i++) {
            frmz[i] = 2*frmz[i] - Pr; 
	}
    }

    while (ctr33 < n) {
	
	ctr22 = start_index;
        while (ctr22 < n_item1) {
            Key t = key(item_array1[ctr22]);
	    if (int_less(t, frmz[ctr33] - tol)) {
	        ctr22++;
                start_index = ctr22; // start_index now points to next item_arry1 element	    
	    }

	    else if (int_eq(frmz[ctr33] - tol, frmz[ctr33] + tol, t)) {
                match = 1;

		// debug
		//printf("Found match %lf\nctr22 = %d\nctr33 = %d\nstart_index = %d\n", *t, ctr22, ctr33, start_index);

                for (int i = 0; i < item_array1[ctr22]->flex_count; i++) {
                    n2_array[item_array1[ctr22]->comb_ind[i]]++;
                    n3_array[item_array1[ctr22]->comb_ind[i]]++;
                }
                *n_item += 1;

                Item *tmp_array = realloc(*item_array, *n_item * sizeof(Item));
                if (!tmp_array) { // Should consider calling a custom exit function to avoid memory leaks on exit (22.06.2020).
                    printf("Error reallocating memory for item_array!\n");
                    free(*item_array);
                    exit(1);
                }
                *item_array = tmp_array;
                (*item_array)[*n_item - 1] = newItem(strlen(item_array1[ctr22]->seq), item_array1[ctr22]->seq, 
				item_array1[ctr22]->enc, item_array1[ctr22]->key, item_array1[ctr22]->fragType, 0);
                sprintf(str4, "%d", (*item_array)[*n_item - 1]->index);
                strcat((*item_array)[*n_item - 1]->fragType, str4);
                if (dblchrg) {
                    strcat((*item_array)[*n_item - 1]->fragType, "(2+)");
                    // need to check if singly charged exists in frmz and if so increment total but this tomorrow - 7.10.22
                    // 03.01.2024: No, this is a mistake, the singly charged fragments will be checked anyway, should remove this.
		    // 17.01.2024: No, it is not a mistake. We are getting n2 > n4 if omitted but will fix it later when back from Tenerife
		    double m1 = 2*frmz[ctr33] - Pr;
                    int index_1 = ctr33 + 1;
                    while (frmz[index_1] <= m1 + tol) { // doesnn't this need ot be <=?
                        if (frmz[index_1] == 0) break;
                        if (int_eq(m1 - tol, m1 + tol, &(frmz[index_1]))) {
                            for (int i = 0; i < item_array1[ctr22]->flex_count; i++) {
                                 n4_array[item_array1[ctr22]->comb_ind[i]]++;
                                 n5_array[item_array1[ctr22]->comb_ind[i]]++;

                            }
                            break;  
                        }
                        index_1++;
                    }
                }
                (*item_array)[*n_item - 1]->intensity = frint[ctr33];
                if (!dblchrg) {
                    
		    int match_h2o = 0;
                    match_h2o = h2o_loss(frmz, ctr33, item_array1[ctr22]->seq, tol, *(item_array1[ctr22]->key), 
				    frint, n5_array, item_array1[ctr22]);
                    
		    if (match_h2o) {
                        //(*n2)++;
                        for (int i = 0; i < item_array1[ctr22]->flex_count; i++) {
                            n3_array[item_array1[ctr22]->comb_ind[i]]++;
                        }
                        *n_item += 1;
                        Item *tmp_array = realloc(*item_array, *n_item * sizeof(Item));
                        if (!tmp_array) { // Should consider calling a custom exit function to avoid memory leaks on exit (22.06.2020).
                            printf("Error reallocating memory for item_array!\n");
                            free(*item_array);
                            exit(1);
                        }
                        *item_array = tmp_array;
                        (*item_array)[*n_item - 1] = newItem(strlen(item_array1[ctr22]->seq), item_array1[ctr22]->seq,
                                item_array1[ctr22]->enc, item_array1[ctr22]->key, item_array1[ctr22]->fragType, 0);

                        *((*item_array)[*n_item - 1]->key) -= H2O;
                        sprintf(str4, "%d", (*item_array)[*n_item - 1]->index);
                        strcat((*item_array)[*n_item - 1]->fragType, str4);
                        strcat((*item_array)[*n_item - 1]->fragType, "-H2O");
                        (*item_array)[*n_item - 1]->intensity = frint[match_h2o];
                    }

                    int match_nh3 = 0;
                    match_nh3 = nh3_loss(frmz, ctr33, item_array1[ctr22]->seq, tol, *(item_array1[ctr22]->key), 
				    frint, n5_array, item_array1[ctr22]);
                    if (match_nh3) {
                        //(*n2)++;
                        for (int i = 0; i < item_array1[ctr22]->flex_count; i++) {
                            n3_array[item_array1[ctr22]->comb_ind[i]]++;
                        }
                        //(*n4)++;
                        // create a new item in item_array
                        *n_item += 1;
                        Item *tmp_array = realloc(*item_array, *n_item * sizeof(Item));
                        if (!tmp_array) { // Should consider calling a custom exit function to avoid memory leaks on exit (22.06.2020).
                            printf("Error reallocating memory for item_array!\n");
                            free(*item_array);
                            exit(1);
                        }
                        *item_array = tmp_array;

                         // now create new item
                        (*item_array)[*n_item - 1] = newItem(strlen(item_array1[ctr22]->seq), item_array1[ctr22]->seq,
                                item_array1[ctr22]->enc, item_array1[ctr22]->key, item_array1[ctr22]->fragType, 0);
			*((*item_array)[*n_item - 1]->key) -= (Nitrogen + 3*H);
                        sprintf(str4, "%d", (*item_array)[*n_item - 1]->index);
                        strcat((*item_array)[*n_item - 1]->fragType, str4);
                        strcat((*item_array)[*n_item - 1]->fragType, "-NH3");
                        (*item_array)[*n_item - 1]->intensity = frint[match_nh3];
                    }
  
                    int match_phos_nloss = 0;
                    match_phos_nloss = phos_nloss(frmz, ctr33, item_array1[ctr22]->seq, tol, 
				    *(item_array1[ctr22]->key), frint, 
				    (*item_array)[*n_item - 1]->enc, n5_array, numMod, item_array1[ctr22]);
                    if (match_phos_nloss) {
                        //(*n2)++;
                        for (int i = 0; i < item_array1[ctr22]->flex_count; i++) {
                            n3_array[item_array1[ctr22]->comb_ind[i]]++;
                        }
                        *n_item += 1;
                        Item *tmp_array = realloc(*item_array, *n_item * sizeof(Item));
                        if (!tmp_array) { // Should consider calling a custom exit function to avoid memory leaks on exit (22.06.2020).
                            printf("Error reallocating memory for item_array!\n");
                            free(*item_array);
                            exit(1);
                        }
                        *item_array = tmp_array;

                        // now create new item
                        (*item_array)[*n_item - 1] = newItem(strlen(item_array1[ctr22]->seq), item_array1[ctr22]->seq,
                                item_array1[ctr22]->enc, item_array1[ctr22]->key, item_array1[ctr22]->fragType, 0);
			*((*item_array)[*n_item - 1]->key) -= 97.99527;
                        sprintf(str4, "%d", (*item_array)[*n_item - 1]->index);
                        strcat((*item_array)[*n_item - 1]->fragType, str4);
                        strcat((*item_array)[*n_item - 1]->fragType, "*");
                        (*item_array)[*n_item - 1]->intensity = frint[match_phos_nloss];
                    }
	        }
		ctr22++;
            }
	    else {
		    
                break;
	    }
	}

	// burned out, will continue tomorrow 02.01.2024

	ctr33++;
    } 

    return match;
    
}


Item int_STsearch(double low, double high) {
    return int_searchR(head, low, high);
}

// Recursive insert.
link insertR(link h, Item item) {
    Key v = key(item);
    //Key t = key(h->item);
    if (h == z) return NEW(item, z, z, 1);
    Key t = key(h->item);
    //if (eq(v, t)) return head;

    if (less(v, t)) h->l = insertR(h->l, item);
    else h->r = insertR(h->r, item);
    (h->N)++;
    return h;  
}

// STinsert midware.
void STinsert(Item item) {
    head = insertR(head, item);
}

// Destroys the tree to free memory.
void STdestroy(link h) {
    // 14.12.2023: possible problem : h should be freed before return
    // let's test. No, it is not that, creates double free() error
    if (h == z) {
        //free(h);
	return;
    }
    //STdestroy(h->r);
    STdestroy(h->l);
    free(h->item->key);
    //free(h->item->comb_ind);
    free(h->item);
    STdestroy(h->r);
    free(h);
    //free(z); // update 15.08.2022	
}

// In order traversal and printing.
void STprint(link h) {
    if (h == z) return;
    STprint(h->l);
    printf("Seq is %s\tmass is %lf\tencoding is %s\tindex is %d\n", h->item->seq, *(h->item->key), h->item->enc, h->item->index);
    STprint(h->r);
   
}

// 31.12.2023: in order traversal and outputing to item array
// item is a pointer to Item struct 
void STto_array(link h, Item *arr, int *i) {
    if (h == z) return;
    STto_array(h->l, arr, i);

    // debug
    //printf("Seq is %s\tmass is %lf\tencoding is %s\tindex is %d\n", h->item->seq, 
    //		    *(h->item->key), h->item->enc, h->item->index);
    
    arr[*i] = h->item;
    (*i)++;
    STto_array(h->r, arr, i);
}

// Range selection function. Uses modified in order traversal.
int rangeSelect(link h, double n1, double n2) {

    // For debugging
    //printf("Range selecting for %lf %lf\n", n1, n2);
    
    if (h == z) return 0;

    // For debugging
    //printf("Seq of head is %s\tmass is %lf\tencoding is %s\ttype is %c\tindex is %d\n", h->item->seq, *(h->item->key), h->item->enc, h->item->fragType, h->item->index);

    /*
    if (*(key(h->item)) >= n1 && *(key(h->item)) <= n2) {
        printf("Seq is %s\tmass is %lf\tencoding is %s\tindex is %d\n", h->item->seq, *(h->item->key), h->item->enc, h->item->index);
    }
    */
    //if (*(key(h->item)) > n1) rangeSelect(h->l, n1, n2); // This does not seem right. Update 15.07.2020: it is right, that how in order traversal works: if n1 > *(key(h->item)) we recur with the left branch. Then we test if it is in range and print it and then recur with right branch.

    if (*(key(h->item)) >= n1 && *(key(h->item)) <= n2) {
        //printf("Search mass is %lf Seq is %s\tmass is %lf\tencoding is %s\ttype is %c\tindex is %d\n", (n2+n1)/2, h->item->seq, *(h->item->key), h->item->enc, h->item->fragType, h->item->index);
        return 1;
    }
    else if (*(key(h->item)) > n2) return rangeSelect(h->l, n1, n2); // NB: need to recur with return. Otherwise always returns 0.
    else if (*(key(h->item)) < n1) return rangeSelect(h->r, n1, n2);
    return 0;
}

/*
// Non-recursive STinsert
void STinsert(Item item) {
    Key v = key(item);
    link p = head, x = p;
    if (head == NULL) {
        head = NEW(item, NULL, NULL, 1);
        return;
    }
    while (x != NULL) {
        p = x;
        x->N++;
        x = less(v, key(x->item)) ? x->l : x->r;
    }
    x = NEW(item, NULL, NULL, 1);
    if (less(v, key(p->item))) p->l = x;
    else p->r = x;   
}
*/

/*
  Now need a function to generate balanced BST from y fragment masses for each generated combination.
  Should use Item type.

  Update 18.07.2020: the tree does not seem balanced and when adding b frags to y frags and frags 
  from different isoforms of the same modified peptide. There are also problems with fragments that 
  have different sets of modifications but the same mass. Need to right 
  proper balanced tree implementation functions. 

  Update 19.07.2020: kinda false alarm yesterday. That is why one should write comments about everything.
  Since I am starting from sorted array of y fragments the tree will grow more or less balanced because 
  all future fragments come from the same sequence but different modifications positions.

  The key is that I use the left<head<=right rule, which handles duplicate keys just fine. Could probably 
  improve performance by implementing self-balancing tree, red-black or avl, but this will require some
  fair amount of extra coding for not so great performance gain.  

  The only thing needed correcting was to account for the fragments that have same mass (key) and sequence
  but different modification sets, which should be rare but anyway. This slows down execution 2-3 fold but
  nothing I can do for the time being.

  Balanced tree algorithm:
  [1,2,3,4,5,6,7,8,9,10]
  1. Get the midlle element and insert it: 10/2 = 5 so it is array[5] = 6
  2. recursively inser left half (array 0 to 5)
  3. The same with right half

                           6
                     .         .     
                    3            9
                   . .          .    .
                  1   4        7       10
                   .   .        . 
                    2   5         8
                                 
                               
                              

   seq is AALSENEQLKK

algorithm start: 
                                  start = 0; end = 10
                                  x = (0 + 10)/2 = 5
                                  item1->seq = NEQLKK

                           .                                           .
                
                 start = 6; end = 10;                            start = 0; end = 4;
                 x = 8;                                          x = (0 + 4)/2 = 2
                 item1->seq = LKK                                item1->seq = LSENEQLKK  

           .                 .                        .                           .

       start 9       start = 6; end = 7;         start = 3; end = 4;               start = 0; end = 1;
       end 10
           x 9       x = 6;                      x = 3;                            x = (0 + 1)/2 = 0
           KK        item1->seq = EQLKK          item1->seq = SENEQLKK             item1->seq =  AALSENEQLKK here is wrong
                                                            
                        
           .        .             .                        .              .                                      .                        .
         
           10      return   start = 7       start = 4;        start = 3; end = 2;         start = 1; end = 1;       start = 0; end = -1; return 
           10               end = 7               end = 4;          return                      x = 1
           K                x = 7;                x = 4;                                        item1->seq = ALSENEQLKK
                          item1->seq = QLKK     item1->seq = ENEQLKK                                   

                                                 .                                     .                       .
                            
                                         return   return                       start = 2; end = 1; return        start = 0; end = 0; return
              

    one extra because it considers precursor as y0

*/

// Helper function that builds BST from sorted data, here with the y fragments of a peptide.
void bst_y(link head, Item item, int start, int end, Item item1) {
    if (start > end || end == 0) return;
    int x = (start + end)/2;

    // Update 15.08.2020: if x == 0 return
    if (x == 0) {
        //bst_y(head, item, start, x - 1, item1);
        //bst_y(head, item, x + 1, end, item1);
        x = 1; // important to take care for start = 0 and end = 1 cases
        //return;
    } 
    
    double extraMass = 0;
    strcpy(item1->seq, item->seq + x);
    strcpy(item1->enc, item->enc + x);
    int str_len = strlen(item1->seq);
    for (int j = 0; j < str_len; j++) {
        if (item1->enc[j] != '.') extraMass = extraMass + ptm[item1->enc[j] - 48].delta_mass;
    }
    strcpy(item1->fragType, "y");
    double Mass = pepMass(item1->seq) + extraMass + Pr; // Important to add proton mass.
    item1->key = &Mass;
    item1->index = x;
    //item1->index = strlen(item1->seq); for y fragments index should be the string length

    // 30.10.2022: need to refactor here. Should capture the return in a variable
    // Item tmp_item = STsearch_ins...
    // if (tmp_item)...
    Item tmp_item = STsearch_ins(item1->key, item1->seq, item1->enc);
    //if (STsearch_ins(item1->key, item1->seq, item1->enc) == NULLitem){
    if (tmp_item == NULLitem) {
        STinsert(newItem(str_len, item1->seq, item1->enc, item1->key, item1->fragType, comb_index));
     
        // for debugging
        //printf("inserted %s\n%s\n", item1->seq, item1->enc);
        
        // Update 21.07.2021: include water and amonia lodsses
        /* Does not seem to improve results but slows down significantly
        if (strstr(item1->seq, "D") || strstr(item1->seq, "E") || strstr(item1->seq, "S") || strstr(item1->seq, "T")){
            double m1 = *(item1->key) - H2O;
            if (STsearch_ins(&m1, item1->seq, item1->enc) == NULLitem){
                STinsert(newItem(x, item1->seq, item1->enc, &m1, item1->fragType, 1, 0));
            }
        }
        if (strstr(item1->seq, "N") || strstr(item1->seq, "Q") || strstr(item1->seq, "K") || strstr(item1->seq, "R")){
            double m1 = *(item1->key) - (Nitrogen + 3*H);
            if (STsearch_ins(&m1, item1->seq, item1->enc) == NULLitem){
                STinsert(newItem(x, item1->seq, item1->enc, &m1, item1->fragType, 0, 1));
            }
        }*/
        
    }
    else {
	// should not grow it before checking if max_flex == flex_count
        /*if (tmp_item->max_flex == tmp_item->flex_count) {
	    Item temp = NULL;
	    //temp = grow_array(tmp_item);
	    temp = realloc(tmp_item, sizeof(*temp) + sizeof(int)*(tmp_item->max_flex + 1));
	    tmp_item = temp;
	}*/
	tmp_item->comb_ind[tmp_item->flex_count++] = comb_index;
	// for debuging
        //printf("fragment exists on BST. encoding is %s, comb_index is %d\n flex_count is %d\n", tmp_item->enc, comb_index, tmp_item->flex_count);
        //for (int i = 0; i < tmp_item->max_flex; i++) printf("comb_ind[%d] = %d\n", i, tmp_item->comb_ind[i]);

    }
    bst_y(head, item, start, x - 1, item1);
    bst_y(head, item, x + 1, end, item1); 
}

// Adding b fragments to the tree...
/*

                                                   seq is AALSENEQLKK
                                                   x = 5
                                                   AALSE 

                             start 0                                       start 6
                             end 4                                         end 10 
                             x 2                                           x 8
                             AA                                            AALSENEQ 

              start 0             start 3                     start 6                    start 9   
              end 1               end 4                       end 7                      end 10
              x 0                 x 3                         x 6                        x 9
              return              AAL                         AALSEN                     AALSENEQL
              here we 
              miss 1                                   start 6     start 7         start 9       start 10
              so x = 1   start 3       start 4         end 5       end 7           end 8         end 10
              A          end 2         end 4           return      x 7             return        x 10
                         return        x 4                         AALSENE                       AALSENEQLK 
         start 0   start 1             AALS
         end 0     end 0
         return    return

*/

void bst_b(link head, Item item, int start, int end, Item item1) {
    if (start > end || end == 0) return;
    int x = (start + end)/2;

    // Update 14.08.2020: need to return if x == 0. This will avoid generating b fragments with index 0.
    if (x == 0){
         x = 1;
    }

    double extraMass = 0;
    strcpy(item1->seq, item->seq);
    strcpy(item1->enc, item->enc);
    item1->seq[x] = '\0';
    item1->enc[x] = '\0';
    strcpy(item1->fragType, "b");

    for (int j = 0; j < x/*strlen(item1->seq)*/; j++) { // Update 15.08.2020: instead of strlen() should use x
        if (item1->enc[j] != '.') extraMass = extraMass + ptm[item1->enc[j] - 48].delta_mass;
    }

    // Update: 07.08.2020: there was an error in b-series mass calculation, need to subtract water
    //double Mass = pepMass(item1->seq) + extraMass + Pr; // Important to add proton as it should be charged.
    double Mass = pepMass(item1->seq) - H2O + extraMass + Pr;
    item1->key = &Mass;
    item1->index = x;
    Item tmp_item = STsearch_ins(item1->key, item1->seq, item1->enc); 
    //if (STsearch_ins(item1->key, item1->seq, item1->enc) == NULLitem){
    if (tmp_item == NULLitem) {
        STinsert(newItem(x, item1->seq, item1->enc, item1->key, item1->fragType, comb_index)); 
        // Update 21.07.2021: include water and amonia losses
        /*
        if (strstr(item1->seq, "D") || strstr(item1->seq, "E") || strstr(item1->seq, "S") || strstr(item1->seq, "T")){
            double m1 = *(item1->key) - H2O;
            if (STsearch_ins(&m1, item1->seq, item1->enc) == NULLitem){
                STinsert(newItem(x, item1->seq, item1->enc, &m1, item1->fragType, 1, 0));
            }
        }
        if (strstr(item1->seq, "N") || strstr(item1->seq, "Q") || strstr(item1->seq, "K") || strstr(item1->seq, "R")){
            double m1 = *(item1->key) - (Nitrogen + 3*H);
            if (STsearch_ins(&m1, item1->seq, item1->enc) == NULLitem){
                STinsert(newItem(x, item1->seq, item1->enc, &m1, item1->fragType, 0, 1));
            }
        }
        */
	// for debugging
        //printf("inserted %s\n%s\n", item1->seq, item1->enc);
    }
    else {
        //if (tmp_item->max_flex == tmp_item->flex_count) tmp_item = grow_array(tmp_item);
	tmp_item->comb_ind[tmp_item->flex_count++] = comb_index;
	// for debuging
        //printf("fragment exists on BST. encoding is %s, comb_index is %d\n flex_count is %d\n", tmp_item->enc, comb_index, tmp_item->flex_count);
        //for (int i = 0; i < tmp_item->max_flex; i++) printf("comb_ind[%d] = %d\n", i, tmp_item->comb_ind[i]);

    
    }
    bst_b(head, item, start, x - 1, item1);
    bst_b(head, item, x + 1, end, item1);
}

// 13.09.2022: there was an error here: should consider theoretical mass, not the detected peak
int nh3_loss(double *frmz, int index, char *seq, double tol, double mass, double *frint, int *n4, Item h) {
    if (index == 0) return 0;
    double m1 = 0;
    if (strstr(seq, "N") || strstr(seq, "Q") || strstr(seq, "K") || strstr(seq, "R")){
        //(*n4)++;
	/*for (int i = 0; i < h->flex_count; i++) {
            n4[h->comb_ind[i]]++;
        }*/
        m1 = mass - (Nitrogen + 3*H);
        while (index > 0) {	
            if (int_eq(m1 - tol, m1 + tol, &(frmz[index]))) {
	        //printf("%lf\n", frmz[index]);
		for (int i = 0; i < h->flex_count; i++) {
                    n4[h->comb_ind[i]]++;
                }
	        return index;
	    }
	    index--;
	}
    }
    return 0;
}

int h2o_loss(double *frmz, int index, char *seq, double tol, double mass, double *frint, int *n4, Item h) {
    if (index == 0) return 0;
    double m1 = 0;
    if (strstr(seq, "D") || strstr(seq, "E") || strstr(seq, "S") || strstr(seq, "T")){
        //(*n4)++;
	/*for (int i = 0; i < h->flex_count; i++) {
            n4[h->comb_ind[i]]++;
        }*/
        m1 = mass - H2O;
        while (index > 0) {
            if (int_eq(m1 - tol, m1 + tol, &(frmz[index]))) {
	        for (int i = 0; i < h->flex_count; i++) {
                    n4[h->comb_ind[i]]++;
                }
                return index;
	    }
            index--;
        }
    }
    return 0;
}

// Phosphorylation neutral loss: to be implemented later (29.09.2022)
int phos_nloss(double *frmz, int index, char *seq, double tol, double mass, double *frint, char *enc, int *n4, int num_Mod, Item h) {
    if (index == 0) return 0;
    double m1 = 0;
    int is_phs = 0;
    int phs_val = 0;

    // there was a bug here, it is not maxMod that counts it is the number of different ptms, which is 
    // determined in main(). The easiest way to refactor is to make this a global extern variable
    // but no, will determine the size of ptm and take it from there
    //long unsigned int numMod = sizeof(ptm)/sizeof(PTM);
   
    // for debugging 
    // 6.10.2022: sizeof would not work because ptm is declared as PTM *ptm, not as PTM ptm[n]. 
    // Will have to refactor this. interval_search should receive the number of different modifications
    // there are in ptm.
    //printf("numMod = %lu\n", numMod);

    for (int i = 0; i < num_Mod; i++) {  // should replace 3 with maxMod but maxMod need to be global
        if (strstr(ptm[i].name, "Phospho")) {
	    is_phs = 1;
            phs_val = i;
	}
    }
    char str1[16];
    sprintf(str1, "%d", phs_val);
    if (is_phs && strstr(enc, str1)){
        m1 = mass - 97.99527;
	//(*n4)++;
	/*for (int i = 0; i < h->flex_count; i++) {
            n4[h->comb_ind[i]]++;
        }*/
        while (index > 0) {
            if (int_eq(m1 - tol, m1 + tol, &(frmz[index]))) {
	        for (int i = 0; i < h->flex_count; i++) {
                    n4[h->comb_ind[i]]++;
                }
	        return index;
	    }
            index--;
        }
    }
    return 0;
}

/* 12.07.2020: How to implement localisation scorring? 
   1. Need to modify BST search to visit all matching fragments and save information in an array. 
      Then use the array to compute scores. The information is in the Item struct. When a fragment
      is matched we need somehow to collect the modification information: there is a loop that steps
      over the enc string of the Item and computes extraMass. We can include there a statement to 
      update a counts[] array for each type of modification. So for each matched fragment we will have
      counts[]. This means we need to alloc and realloc an array of counts[] data structure. Then combine
      any 2 of the counts[] and compute bayesian probabilities. Could it be done in a simpler whay?

*/ 
