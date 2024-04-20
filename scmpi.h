#ifndef scmpi_h__
#define scmpi_h__

#include "sqlite3.h"

#define maxSize 210

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

// Parser function for parsing run parameters
extern int parse_params(char *config_file, int *prTol, double *frTol, char fasta[5][2048], int *maxMod, int *n_losses);
extern double get_Pscore (int, int, double);
extern int scorer(spectrum, results*, results*, int, sqlite3*, Item,
                Item, PTM*, double, double, double [][100][140]);
extern void make_score_table(int n2, int n4, int q2, double table[][n2][n4]);

extern int seq_count;
extern int maxCount;
extern sqlite3 *db;
extern PTM *ptm;
extern int precTol;
extern double tol;
FILE *res;


#endif  // scmpi_h__
