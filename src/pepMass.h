#ifndef PEPMASS_H
#define PEPMASS_H

// Define chemical masses (from ionsource.com).
#define C  12.0
#define Nitrogen 14.003074
#define O  15.99491461986//15.99491463
#define H  1.00782503207 //1.007825035
#define Pr 1.00727646677 //1.00727646688
#define H2O 18.010564684
//#define tol 0.5 // this should be a parameter taken from config file
//#define precTol 5 // as above for tol
#define maxSize 210 // as above for tol
#define q 10.0 // as above for tol; it is important to be 5.0, otherwise (if 5) scorring does not work

//extern int precTol;
//extern double tol;
//extern int q;
//extern int maxSize;


// Struct to hold spectrum data.
typedef struct Spectrum {
    int scan;
    int scanNum;
    int charge;
    double retTime;
    double precMass;
    double fragMass[maxSize];
    double fragInt[maxSize];
    int actual_size;
} spectrum;

// Struct to hold modifications data.
typedef struct { 
    char name[1024];
    double delta_mass;
    char site[16];
    int count;
} PTM;
 
// For debugging. Update 18.07.2020: not for debugging, this is how it should be so all functions see it
// update 19.08.2021: changing it now so openMP can work
extern PTM *ptm;
extern int maxMod;
#pragma omp threadprivate(ptm)

// Struct to hold modification indices. 
typedef struct mod_Index {
    int index;
    int type;
} mod_ind;

static const double aa_table[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 71.03711, 0, 103.00919, 115.02694, 129.04259,
            147.06841, 57.02146, 137.05891, 113.08406, 0, 128.09496,
            113.08406, 131.04049, 114.04293, 114.07931, 97.05276, 128.05858, 156.10111,
            87.03203, 101.04768, 0, 99.06841, 186.07931, 0, 163.06333};

// A function to return AA mass.
double AminoAcidMass(char);

// A function to calculate peptide mass.
double pepMass(char*);

#endif
