#ifndef ITEM_OPEN_H
#define ITEM_OPEN_H

// Struct to hold fragments. 
typedef double *Key;
typedef struct {
    int index;
    char seq[1024];
    char enc[1024];
    Key key;
    char fragType[1024];
    int matched; // 16.09.2022: major refactoring of magvb to make scoring better
    int max_flex;
    int flex_count;
    double intensity; // 29.09.2022: need this to be able to do intensity scoring
    int comb_ind[]; // 30.10.2022 flexible array to keep track of isoform matches
} *Item;

//typedef struct Frag *Item;

#define NULLitem NULL
/*
// Struct to hold the site probabilities 
typedef struct {
    int mod_type;
    int prob[1024];
} site_prob;
*/
// Functions to test keys relations.
int eq(Key, Key);
int less(Key, Key);
int int_eq(double, double, Key);
int int_less(Key, double);
// A function to extract key from item.
Key key(Item);

// A function to create new item.
Item newItem(int, char*, char*, Key, char*, int);

// Update 18.07.2020: new function for comparing items based on key, seq and modifications
int item_comp(Item, char*);
int cmpfunc( const void *, const void *);
#endif
