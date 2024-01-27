/*
 * 1. Reads first search results 
 * 2. Fits model
 * 3. Returns array of coefficients
 * need to parse the result file to 3 arrays
 * PrecMass, RetTime, MassError
 */

#include <stdio.h>
#include <gsl/gsl_multifit.h>
#include <string.h>

int mass_recalibration(int i, char *f, double *res, double min_score, int max_data) {

    /*
     * arguments :
     * 1. number of different PTM considered
     * 2. data file name
     * 3. Array name to save results 
     * 4. Minimal score 
     * 5. Max number of data points to use
     */

    // 2d array to hold tokens during parcing
    char tokens[16 + 2*i][2048];
    char data_str[4096];
    
    double PrecMass_arr[max_data];
    double RetTime_arr[max_data];
    double MassError_arr[max_data];
    
    FILE *d_fl;
    // open data file for reading
    if ((d_fl = fopen(f, "r")) == NULL) {
        printf("Error opening data file!\n");
        exit(1);
    }

    int count_arrays = 0;
    // Tokenize file lines
    while (fgets(data_str, sizeof(data_str), d_fl) != NULL) {
        char *token = strtok(data_str, "\t");
        int t = 0;
        while (token != NULL) {
            strcpy(tokens[t++], token);
            token = strtok(NULL, "\t");
        }

	// now line is stored in tokens row and we can extract the arrays values
        if (atof(tokens[15]) >= min_score) {
            // will fill the arrays until last available slot
            if (count_arrays < max_data) {
                PrecMass_arr[count_arrays] = atof(tokens[13]);
	        RetTime_arr[count_arrays] = atof(tokens[10]);
	        MassError_arr[count_arrays] = atof(tokens[14]);
	        count_arrays++;
	    }
	    // then will replace last elements if there are more elements to fill
	    else {
                PrecMass_arr[count_arrays - 1] = atof(tokens[13]);
                RetTime_arr[count_arrays - 1] = atof(tokens[10]);
                MassError_arr[count_arrays - 1] = atof(tokens[14]);
	    }
	}
    }

    // debug
    printf("Collected %d datapoints with Score > %lf\n", count_arrays, min_score);
    
    // now fit model
    int j, n;
    n = count_arrays;
    double chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;
    
    X = gsl_matrix_alloc (n, 4);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);
    
    c = gsl_vector_alloc (4);
    cov = gsl_matrix_alloc (4, 4);

    for (j = 0; j < n; j++) {
        gsl_matrix_set (X, j, 0, 1.0);
        gsl_matrix_set (X, j, 1, PrecMass_arr[j]);
        gsl_matrix_set (X, j, 2, RetTime_arr[j]);
        gsl_matrix_set (X, j, 3, PrecMass_arr[j]*RetTime_arr[j]);

        gsl_vector_set (y, j, MassError_arr[j]);
    }

    {
        gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 4);
        gsl_multifit_linear (X, y, c, cov, &chisq, work);
        gsl_multifit_linear_free (work);
    }

#define Cc(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  {
    res[0] = Cc(0);
    res[1] = Cc(1);
    res[2] = Cc(2);
    res[3] = Cc(3);

    printf ("# best fit: Y = %g + %g PrecMass + %g RetTime + %g PrecMass*RetTime\n",
            Cc(0), Cc(1), Cc(2), Cc(3));

    printf ("# covariance matrix:\n");
    printf ("[ %+.5e, %+.5e, %+.5e, %+.5e  \n",
               COV(0,0), COV(0,1), COV(0,2), COV(0,3));
    printf ("  %+.5e, %+.5e, %+.5e, %+.5e  \n",
               COV(1,0), COV(1,1), COV(1,2), COV(1,3));
    printf ("  %+.5e, %+.5e, %+.5e, %+.5e  \n",
               COV(2,0), COV(2,1), COV(2,2), COV(2,3));
     printf ("  %+.5e, %+.5e, %+.5e, %+.5e ]\n",
               COV(3,0), COV(3,1), COV(3,2), COV(3,3));
    printf ("# chisq = %g\n", chisq);
  }

  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (w);
  gsl_vector_free (c);
  gsl_matrix_free (cov);

  return 0;

}
