/*
  A short library that parses MS2 structs with variable q 
  Implements 100Th interval 
  limit of the allowed number of peaks.
  Usage: functions are called from callback in scorer_omp
                            
  Metodi V. Metodiev, 2022, Colchester, UK 
*/

#include <stdio.h>
#include <stdlib.h> // for the exit() function
#include <string.h>
#include <ctype.h> // for isdigit() function
#include <math.h>
#include "pepMass.h"
#include "parse_spectrum.h"

int intervalFill(double mass, double intensity, spectrum *spec, int *i, int interval, int q1) {
    double minInt = 0;
    int sorted = 0; // Need to change it somehow to avoid always sorting the array. Either make sorted global or return its value.
 
    if (i[interval] < q1) {
        spec->fragMass[(interval - 1)*q1 + i[interval]] = mass;
        spec->fragInt [(interval - 1)*q1 + i[interval]] = intensity;
        //printf("%lf\n", spec->fragMass[(interval - 1)*q + i[interval]]);
        return 0;
    }
    else {
        
        int c = 0;
        int d = 0;
        double p = 0;
        double r = 0;

        // Sort the interval now if not sorted. Needs to be fixed.
        if (sorted == 0){
            for (c = (interval - 1)*q1; c < interval*q1; c++) {
                d = c;
                while (d > (interval - 1)*q1 && spec->fragInt[d-1] > spec->fragInt[d]) {
                    p = spec->fragMass[d];
                    r = spec->fragInt[d];
                    spec->fragMass[d] = spec->fragMass[d - 1];
                    spec->fragInt[d] = spec->fragInt[d - 1];
                    spec->fragMass[d - 1] = p;
                    spec->fragInt[d - 1] = r;
                    d--;
                }
            }
            sorted = 1;
            // Assign first element to minInt.
            minInt = spec->fragInt[(interval - 1)*q1];
            // If new fragments have higher intensity replace first element with it.
            if (intensity > minInt) {
                spec->fragMass[(interval - 1)*q1] = mass;
                spec->fragInt[(interval - 1)*q1] = intensity;
                sorted = 0;
                return 0;
            }
        }
    } 

    return 0;      
} 

// final sort to make sure fragments are sorted in increasing mass order
void sort_masses(spectrum *spec) {
    int c = 0;
    int d = 0;
    double p = 0;
    double r = 0;
    for (c = 0; c < spec->actual_size; c++) {
        d = c;
	while (d > 0 && spec->fragMass[d-1] > spec->fragMass[d]) {
            p = spec->fragMass[d];
            r = spec->fragInt[d];
            spec->fragMass[d] = spec->fragMass[d - 1];
            spec->fragInt[d] = spec->fragInt[d - 1];
            spec->fragMass[d - 1] = p;
            spec->fragInt[d - 1] = r;
            d--;
       }



    }
}

void remove_0(spectrum *spec) {
    // now remove zeroes
    double tmp_mass[maxSize];
    double tmp_int[maxSize];
    for (int i = 0; i < maxSize; i++) {
        tmp_mass[i] = spec->fragMass[i];
        tmp_int[i] = spec->fragInt[i];
    }
    for (int i = 0; i < maxSize; i++) {
        spec->fragMass[i] = 0.0;
        spec->fragInt[i] = 0.0;
    }
    // now copy only non zero elements to spec.fragMass and spec.fragInt
    int j = 0;
    for (int i = 0; i < maxSize; i++) {
        if (tmp_mass[i] != 0) {
            spec->fragMass[j] = tmp_mass[i];
            spec->fragInt[j] = tmp_int[i];
            j++;
        }
    }
    spec->actual_size = j;
    sort_masses(spec);
}

