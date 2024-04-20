/*
 *
 * Update 23.09.2021: the program runs and debugs but there are problems:
 * 1. Results on mac os and linux cluster differ
 * 2. Large arrays are declared automatic, which could cause stack overflow
 * 4. No way to know that messages are recvd or send before corresponding
 * data is used
 * Should consider:
 * 1. Allocating the spectra and results array on the heap in each process
 * 2. Change to non-blocking communications
 * 3. Write some rest to make sure slave processes are copying the database
 *    properly
 *
 *
 * Update 30.01.2022: there is still the problem with overlapping specificities. 
 * The mod_comb table contains a lot of these. The algorithm does not handle
 * this and the probgram is not functional for now. Should either try to fix
 * it or publish without the open search for combinatorial search first.
 *
 * Update 16.02.2022: tests are promissing but there are certain things
 * that need addresssing:
 * 1. The error tolerance for mod_comb is too large and sulfonation and
 *    phosphorylation are equally selected. Will change it to 10ppm
 * 2. Need to correct scorer_omp to save fragment masses and then 
 *    correct select_by_prob to process the new format. As of now
 *    it does not process socrer_mpi output. In fact it writes huge
 *    files that can fill the storage quickly. This is a bug to carefully
 *    analyse and make sure it is corrected.
 * 3. change the code to determine number of spectra as in scorer_omp.
 *
 *
 * Update 30.10.2022: changing to separate BSTs for isoforms now
 * Will implement flexible array in Item type to keep track of isoforms.
 * Need to refactor code to properly allocate and destroy items and item_array.
 *
 * 20.04.2024: will try to make it portable by writing a wrapper. This is to
 * address problems with newer/older version of ompi at targeted platforms.
*/


#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "sqlite3.h"
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "scmpi.h"

sqlite3 *db;

int main(int argc, char **argv) {

    MPI_Status status;
    spectrum spec1;
    //results result;
    int N_spectra = atoi(argv[4]);


    // mpi variables
    int my_id, num_procs, root_process, avg_spectra_per_process, an_id,
    start_spec, end_spec, num_spec_to_send, num_spec_to_receive;

    // make mpi data types for spectrum and results
    MPI_Datatype spectrum_type;
    MPI_Datatype type1[8] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    int blocklen1[8] = {1, 1, 1, 1, 1, maxSize, maxSize, 1};
    MPI_Aint disp1[8];

    MPI_Datatype results_type;
    MPI_Datatype type2[25] = {MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR,
                              MPI_CHAR, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_CHAR,
                              MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR, MPI_CHAR, MPI_CHAR,
                              MPI_CHAR, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
    int blocklen2[25] = {1, 1, 2048, 2048, 2048, 4096, 4096, 1, 1, 1, 1024, 1, 1, 1, 4096,
                         1024, 1024, 1024, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint disp2[25];


    // process spectra here
    //int ierr = 0;
    MPI_Init(&argc, &argv);
    int ierr = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

     // testing the score lookup table
    static double score_table[6][100][140];
    make_score_table(100, 140, 6, score_table);
    printf("Testing the score lookup table: %lf, %lf, %lf\n",  score_table[0][20][50],
                    score_table[3][40][42], score_table[5][20][22]);

    // allocate arrays
    spectrum *spectra = (spectrum *)calloc(N_spectra, sizeof(spectrum));
    if (!spectra) {
        printf("Error allocating memory for spectra array!\n");
        exit(1);
    }
    results *results1 = (results *)calloc(N_spectra, sizeof(results));
    if (!spectra) {
        printf("Error allocating memory for results1 array!\n");
        exit(1);
    }
    results *results2 = (results *)calloc(N_spectra, sizeof(results));
    if (!spectra) {
        printf("Error allocating memory for results2 array!\n");
        exit(1);
    }

    // now create item, item1, and ptm
    //PTM ptm[3];
    Item item = NULL;
    Item item1 = NULL;

    // commit MPI_datatypes
    disp1[0] = (char *)&spectra[0].scan - (char *)&spectra[0];
    disp1[1] = (char *)&spectra[0].scanNum - (char *)&spectra[0];
    disp1[2] = (char *)&spectra[0].charge - (char *)&spectra[0];
    disp1[3] = (char *)&spectra[0].retTime - (char *)&spectra[0];
    disp1[4] = (char *)&spectra[0].precMass - (char *)&spectra[0];
    disp1[5] = (char *)&spectra[0].fragMass - (char *)&spectra[0];
    disp1[6] = (char *)&spectra[0].fragInt- (char *)&spectra[0];
    disp1[7] = (char *)&spectra[0].actual_size - (char *)&spectra[0];
    MPI_Type_create_struct(8, blocklen1, disp1, type1, &spectrum_type);
    MPI_Type_commit(&spectrum_type);

    disp2[0] = (char *)&results1[0].Total - (char *)&results1[0];
    disp2[1] = (char *)&results1[0].Matched - (char *)&results1[0];
    disp2[2] = (char *)&results1[0].mod1_probs - (char *)&results1[0];
    disp2[3] = (char *)&results1[0].mod2_probs - (char *)&results1[0];
    disp2[4] = (char *)&results1[0].mod3_probs - (char *)&results1[0];
    disp2[5] = (char *)&results1[0].Fragments - (char *)&results1[0];
    disp2[6] = (char *)&results1[0].Masses - (char *)&results1[0];
    disp2[7] = (char *)&results1[0].Scan - (char *)&results1[0];
    disp2[8] = (char *)&results1[0].PrecScan - (char *)&results1[0];
    disp2[9] = (char *)&results1[0].RetTime - (char *)&results1[0];
    disp2[10] = (char *)&results1[0].Sequence - (char *)&results1[0];
    disp2[11] = (char *)&results1[0].Charge - (char *)&results1[0];
    disp2[12] = (char *)&results1[0].PrecMass - (char *)&results1[0];
    disp2[13] = (char *)&results1[0].Score - (char *)&results1[0];
    disp2[14] = (char *)&results1[0].Proteins - (char *)&results1[0];
    disp2[15] = (char *)&results1[0].mod1 - (char *)&results1[0];
    disp2[16] = (char *)&results1[0].mod2 - (char *)&results1[0];
    disp2[17] = (char *)&results1[0].mod3 - (char *)&results1[0];
    disp2[18] = (char *)&results1[0].nMod - (char *)&results1[0];
    disp2[19] = (char *)&results1[0].Missed - (char *)&results1[0];
    disp2[20] = (char *)&results1[0].Decoy - (char *)&results1[0];
    disp2[21] = (char *)&results1[0].Contam - (char *)&results1[0];
    disp2[22] = (char *)&results1[0].seq_count - (char *)&results1[0];
    disp2[23] = (char *)&results1[0].score_threshold - (char *)&results1[0];
    disp2[24] = (char *)&results1[0].matched_int - (char *)&results1[0];

    MPI_Type_create_struct(25, blocklen2, disp2, type2, &results_type);
    MPI_Type_commit(&results_type);

    root_process = 0;

    // allocate memory for item and item1
    // 31.10.2022: changing it to accomodate flexible arrays. Will allocate 100 elements initially.
    // Then at each operation will test if n_flex is < then max_flex and reallocate if necessary.
    item = malloc(sizeof(*item) + 1*sizeof(int));  // This might need to go to scorer_mod.c and be passed as argument.
    if (!item) {
        printf("Error allocating memory for item!\n");
        exit(1);
    }

    // for debugging
    //printf("Created item in thread %d address is %x\n", omp_get_thread_num(), &item);

    item1 = malloc(sizeof(*item1) + 1*sizeof(int)); // This might need to go to scorer_mod.c and be passed as argument.
    if (!item1) {
        printf("Error allocating memory for item!\n");
        exit(1);
    }

    ptm = (PTM *) malloc(3 * sizeof(PTM));
    if (!ptm) {
        printf("Error allocating memory for ptm array!\n");
        exit(1);
    }

    //sqlite3 *db; // needs to be global
    int rc;
    sqlite3_stmt *ppStmt;
    rc = sqlite3_open(argv[2], &db);
    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return(1);
    }

    // copy to in memory database
    // now need to copy some data to an in-memory db that will be used in eval_obj()
    sqlite3 *db1;
    sqlite3_backup *pBackup;


    rc = sqlite3_open(":memory:", &db1); // need to write error handling code later
    pBackup = sqlite3_backup_init(db1, "main", db, "main");
    rc = sqlite3_backup_step(pBackup, -1);
    if (rc == SQLITE_DONE) {
        printf("Finished copying to in-memory database in process %d\n", my_id);
    }
    else {
        fprintf(stderr, "Can't copy to in memory database: %s\n", sqlite3_errmsg(db1));
        sqlite3_close(db1);
        return(1);
    }
    (void)sqlite3_backup_finish(pBackup);


    // update 27.09.2021: implementing load balancing
    // will have to use scatterv and gatherv
    avg_spectra_per_process = N_spectra / num_procs;
    int n_extra = N_spectra % num_procs;
    int k = 0;
    int send_count[num_procs];
    int recv_count = 0;
    int displ[num_procs];
    for (int i = 0; i < num_procs; i++) {
        if (i < n_extra) send_count[i] = avg_spectra_per_process + 1;
        else send_count[i] = avg_spectra_per_process;
        displ[i] = k;
        k = k + send_count[i];
    }
    recv_count = send_count[my_id];

    for (int t = 0; t < N_spectra; t++) {
        results1[t].Score = 0;
        results2[t].Score = 0;
        strcpy(results2[t].Fragments, "");
        strcpy(results1[t].Fragments, "");
        strcpy(results1[t].Masses, "");
        strcpy(results2[t].Masses, "");
        strcpy(results1[t].mod1_probs, "");
        strcpy(results2[t].mod1_probs, "");
        strcpy(results1[t].mod2_probs, "");
        strcpy(results2[t].mod2_probs, "");
        strcpy(results1[t].mod3_probs, "");
        strcpy(results2[t].mod3_probs, "");

    }

    if (my_id == root_process) {


        char fasta[5][2048]; // up to 5 fasta files can be used
        int maxMod = 0;
        int n_losses = 0;
        int ct1 = 0;




        /*for (int t = 0; t < N_spectra; t++) {
            results1[t].Score = 0;
            results2[t].Score = 0;
            strcpy(results2[t].Fragments, "");
            strcpy(results1[t].Fragments, "");
            strcpy(results1[t].Masses, "");
            strcpy(results2[t].Masses, "");
            strcpy(results1[t].mod1_probs, "");
            strcpy(results2[t].mod1_probs, "");
            strcpy(results1[t].mod2_probs, "");
            strcpy(results2[t].mod2_probs, "");
            strcpy(results1[t].mod3_probs, "");
            strcpy(results2[t].mod3_probs, "");

        }*/

        // File pointer to read config info.
        FILE *cnfg;
        if ((cnfg = fopen(argv[3], "r")) == NULL) {
            printf("Error opening config file for reading.\n");
            exit(1);
        }

        // Parse config file
        int fasta_cntr = parse_params(argv[3], &precTol, &tol, fasta, &maxMod, &n_losses);

        // for debugging
        printf("precTol = %d\n", precTol);
        printf("tol = %lf\n", tol);
        printf("maxMod = %d\n", maxMod);
        printf("FASTA files are:\n");
        for (int j = 0; j < fasta_cntr; j++) printf("%s\n", fasta[j]);

        MPI_Bcast(&precTol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


        FILE *mms_ptr;
        // Read data struct by struct, search and write results to file.
        if ((mms_ptr = fopen(argv[1], "rb")) == NULL) {
            printf("Error opening binary file for reading.\n");
            exit(1);
        }
        char resFile[128] = {};
        strcat(resFile, argv[1]);
        strcat(resFile, ".txt");

        if ((res = fopen(resFile, "w")) == NULL) {
            printf("Error opening results file for writing.\n");
            exit(1);
        }

        fprintf(res, "Total\tMatched"); // Need matched and total ion numbers (18.06.2020).
        for (int ctr = 0; ctr < 3; ctr++) fprintf(res, "\tmod%d probs", ctr + 1); // later should put results in data struct perhaps
        fprintf(res, "\tFragments\tMasses\tScan\tPrecScan\tRetTime\tSequence\tCharge\tPrecMass\tScore\tProteins");
        for (int ctr = 0; ctr < 3; ctr++) fprintf(res, "\tmod%d", ctr + 1);
        fprintf(res, "\tnMod\tMissed\tDecoy\tContam\tSeq_count\tscore_threshold\tmatched_int\n");

        // Update 16.02.2022: will count number of spectra and increment N_spectra
        /*while (fread(&spec1, sizeof(spectrum), 1, mms_ptr) == 1) {
            N_spectra++;
        }
        // will not do it for now as it will require major refactoring
        */
        for (int j = 0; j < N_spectra; j++) {
            // this should be only done by master
            fseek(mms_ptr, sizeof(spectrum) * j, SEEK_SET);
            fread(&spectra[j], sizeof(spectrum), 1, mms_ptr);

            /*printf("Prec mass is %lf\n", spectra[j].precMass);
            printf("Charge is %d\n", spectra[j].charge);
            printf("fragMass[10] is %lf\n", spectra[j].fragMass[10]);*/


        }


        // scatterv the spectra to all processes
        MPI_Scatterv(&spectra[0], send_count, displ, spectrum_type, &spectra[0], recv_count, spectrum_type, 0, MPI_COMM_WORLD);

        // for debugging
        /*for (int i = 0; i < recv_count; i++) {
            printf("precMass[%d] at process %d = %lf\n", i, my_id, spectra[i].precMass);
            printf("charge[%d] at process %d = %d\n", i, my_id, spectra[i].charge);
            printf("fragMass[%d] at process %d = %lf\n", i, my_id, spectra[i].fragMass[10]);
        }*/
        printf("Root will process the first %d spectra\n", send_count[0]);

        for (int w = 0; w < send_count[0]; w++) {
            //printf("Prec mass is %lf\n", spectra[w].precMass);
            scorer(spectra[w], results1, results2, w, db1, item, item1, ptm, precTol, tol, score_table);
            // 29.10.2022: adding speq_count
	    results1[w].seq_count = seq_count;
	    results2[w].seq_count = seq_count;
	    /*if (results1[w].Total != 0) results1[w].score_threshold = -10*log10(0.01/seq_count);
	    results2[w].seq_count = seq_count;
	    if (results2[w].Total != 0) results2[w].score_threshold = -10*log10(0.01/seq_count);
	    seq_count = 0;*/
    	    //printf("total at process %d is %d\n", my_id, results1[w].Total);

            // 03.03.2024: seq_count need to be reset
            seq_count = 0;
	}


        // receive data from slaves
        MPI_Gatherv(&results1[0], recv_count, results_type, &results1[0], send_count, displ, results_type, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&results2[0], recv_count, results_type, &results2[0], send_count, displ, results_type, 0, MPI_COMM_WORLD);

        /*for (int w = 0; w < recv_count; w++) {
            printf("total at process %d is %d\n", my_id, results1[w].Total);
        }*/

        // ok, it is working in the master thread, now make it work with slaves

        for (int i = 0; i < N_spectra; i++) {
            if (results1[i].Total != 0) results1[i].score_threshold = -10*log10(0.01/results1[i].seq_count);
            if (results2[i].Total != 0) results2[i].score_threshold = -10*log10(0.01/results2[i].seq_count);

            fprintf(res, "%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%lf\t%s\t%d\t%lf\t%lf\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n",\
        results1[i].Total, results1[i].Matched, results1[i].mod1_probs, results1[i].mod2_probs, \
        results1[i].mod3_probs, results1[i].Fragments, results1[i].Masses, results1[i].Scan, results1[i].PrecScan, \
        results1[i].RetTime, results1[i].Sequence, results1[i].Charge, results1[i].PrecMass, \
        results1[i].Score, results1[i].Proteins, results1[i].mod1, results1[i].mod2, results1[i].mod3, \
        results1[i].nMod, results1[i].Missed, results1[i].Decoy, results1[i].Contam, results1[i].seq_count, results1[i].score_threshold, results1[i].matched_int);

            fprintf(res, "%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%lf\t%s\t%d\t%lf\t%lf\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n",\
        results2[i].Total, results2[i].Matched, results2[i].mod1_probs, results2[i].mod2_probs, \
        results2[i].mod3_probs, results2[i].Fragments, results2[i].Masses, results2[i].Scan, results2[i].PrecScan, \
        results2[i].RetTime, results2[i].Sequence, results2[i].Charge, results2[i].PrecMass, \
        results2[i].Score, results2[i].Proteins, results2[i].mod1, results2[i].mod2, results2[i].mod3, \
        results2[i].nMod, results2[i].Missed, results2[i].Decoy, results2[i].Contam, results2[i].seq_count, results2[i].score_threshold, results2[i].matched_int);

        }

        fclose(mms_ptr);
        fclose(cnfg);
        fclose(res);


    }

    else {

        MPI_Bcast(&precTol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Scatterv(&spectra[0], send_count, displ, spectrum_type, &spectra[0], recv_count, spectrum_type, 0, MPI_COMM_WORLD);

        // for debugging
        /*for (int i = 0; i < recv_count; i++) {
            printf("precMass[%d] at process %d = %lf\n", i, my_id, spectra[i].precMass);
            printf("charge[%d] at process %d = %d\n", i, my_id, spectra[i].charge);
            printf("fragMass[%d] at process %d = %lf\n", i, my_id, spectra[i].fragMass[10]);
        }*/

        // initialize the results arrays
        /*for (int t = 0; t < N_spectra; t++) {
            results1[t].Score = 0;
            results2[t].Score = 0;
            strcpy(results2[t].Fragments, "");
            strcpy(results1[t].Fragments, "");
            strcpy(results1[t].Masses, "");
            strcpy(results2[t].Masses, "");
            strcpy(results1[t].mod1_probs, "");
            strcpy(results2[t].mod1_probs, "");
            strcpy(results1[t].mod2_probs, "");
            strcpy(results2[t].mod2_probs, "");
            strcpy(results1[t].mod3_probs, "");
            strcpy(results2[t].mod3_probs, "");

        }*/


        for (int w = 0; w < recv_count; w++) {
            scorer(spectra[w], results1, results2, w, db1, item, item1, ptm, precTol, tol, score_table);
            //printf("total at process %d is %d\n", my_id, results1[w].Total);
	    results1[w].seq_count = seq_count;
	    results2[w].seq_count = seq_count;
	    seq_count = 0;
        }

        // for debuging
        /*for (int w = 0; w < recv_count; w++) {
            printf("total at process %d is %d\n", my_id, results1[w].Total);
        }*/


        MPI_Gatherv(&results1[0], recv_count, results_type, &results1[0], send_count, displ, results_type, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&results2[0], recv_count, results_type, &results2[0], send_count, displ, results_type, 0, MPI_COMM_WORLD);
    }

    sqlite3_close(db);
    sqlite3_close(db1);

    free(item);
    free(item1);
    free(ptm);
    free(spectra);
    free(results1);
    free(results2);

    // finalize
    ierr = MPI_Finalize();





    //printf("Processed %d MS/MS spectra.\n", ct1);

    return 0;
}

