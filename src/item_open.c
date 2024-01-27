/* 

  Definitions of Item ADT 

*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "item_open.h"
#include "pepMass.h"

int eq(Key a, Key b) {
    if (*a == *b) return 1;
    else return 0;
}

int int_eq(double low, double high, Key key) {
    if (*key >= low && *key <= high) return 1;
    else return 0;
}

int less(Key a, Key b) {
    if (*a < *b) return 1;
    else return 0;   
}

int int_less(Key key, double low) {
    if (*key < low) return 1;
    else return 0;    
}

Key key(Item item) {
    return item->key;
}

Item newItem(int index, char* seq, char* enc, Key key, char *fragType, int comb_index) {
    Item new = malloc(sizeof *new + 100*sizeof(int));
    new->index = index;
    new->key = malloc(sizeof *key); // Why do we need malloc here?
    *(new->key) = *key;
    strcpy(new->seq, seq);
    strcpy(new->enc, enc);
    strcpy(new->fragType, fragType);
    new->matched = 0; // 16.09.2022
    new->max_flex = 100;
    for (int i = 0; i < new->max_flex; i++) new->comb_ind[i] = 0;
    new->flex_count = 1;
    new->comb_ind[0] = comb_index;
    return new;
}

int cmpfunc( const void *a, const void *b) {
    return *(char*)a - *(char*)b;
}

// Compares two items by their keys, sequence and modification sets
/* Update 20.07.2020; It is givin errors: some strings that have the same modification sets 
   are evaluated as different. This causes it to use GB memory.

unsorted item->enc: ......0.........
sorted item->enc:
unsorted enc: ..............0.
sorted enc: ...............0
return = -46
Found it!
New seq is	DSKPSSSSPVAYVETK
New enc is	..............0.
h->item->seq is	DSKPSSSSPVAYVETK
h->item->enc is	......0.........

The item->enc is lost after sorting.

Solved: it turns out strcpy should be given array[1024] instead of arr[strlen()].

01.02.2024: this is a bottleneck function. Should replace it with something like:
look at the ptms in Item abd is they are the same return 0.
*/
int item_comp_new(int n_mod, Item item, int counts[n_mod], int i, PTM *ptm) {
    int res = 1;
    for (int j = 0; j < i; j++) {
        if (counts[j] == ptm[j].count) continue;
        else res = 0;	
    }
    return res;
}


int item_comp(Item item, char *enc) {

/*    int cmpfunc( const void *a, const void *b) {
        return *(char*)a - *(char*)b;
      }
*/
    int size1 = strlen(item->enc);
    int size2 = strlen(enc);
    char item1_enc[1024];
    char item2_enc[1024];
 
    strcpy(item1_enc, item->enc);
    strcpy(item2_enc, enc);
 
    qsort(item1_enc, size1, sizeof(char), cmpfunc);
    qsort(item2_enc, size2, sizeof(char), cmpfunc);

    /* For debugging: it does not find any with maxMod 6 and phos, acetyl and metOx. Should check from time to time 
       with other modifications and if never finds any should probably drop the item_com call in searchR_ins().
    printf("unsorted item->enc: %s\n", item->enc);
    printf("sorted item->enc: %s\n", item1_enc);
    printf("unsorted enc: %s\n", enc);
    printf("sorted enc: %s\n", item2_enc);
    printf("return = %d\n", strcmp(item1_enc, item2_enc)); 
    */
    return strcmp(item1_enc, item2_enc);
}

