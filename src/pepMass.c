#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pepMass.h"

// Definition of ptm declared in the header file as extern.
PTM *ptm = NULL;
int maxMod = 0;
#pragma omp threadprivate(ptm)

//int precTol = 10;
//double tol = 0.6;

// A function to return AA mass.
double AminoAcidMass(char c) {
    // 13.01.2024: will make it into lookup table to see if it speeds up
    /*static const double aa_table[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	    0, 0, 0, 0, 71.03711, 0, 103.00919, 115.02694, 129.04259,
	    147.06841, 57.02146, 137.05891, 113.08406, 0, 128.09496,
	    113.08406, 131.04049, 114.04293, 114.07931, 97.05276, 128.05858, 156.10111,
	    87.03203, 101.04768, 0, 99.06841, 186.07931, 0, 163.06333};
   */
   return aa_table[c - 42];
   /*switch (c) {
      case 'A': return  71.03711;
      case 'C': return 103.00919;
      case 'D': return 115.02694;
      case 'E': return 129.04259;
      case 'F': return 147.06841;
      case 'G': return  57.02146;
      case 'H': return 137.05891;
      case 'I': return 113.08406;
      case 'K': return 128.09496;
      case 'L': return 113.08406;
      case 'M': return 131.04049;
      case 'N': return 114.04293;
      case 'O': return 114.07931;
      case 'P': return  97.05276;
      case 'Q': return 128.05858;
      case 'R': return 156.10111;
      case 'S': return  87.03203;
      case 'T': return 101.04768;
      case 'V': return  99.06841;
      case 'W': return 186.07931;
      case 'Y': return 163.06333;
      case '*': return 0; // might not be necessary though.
    }
  return 0;*/
}

// A function to calculate peptide mass.
double pepMass(char *p) {
    int strl = strlen(p);
    if (strl == 0) return 0;
    //printf("Sequence is %s\n", p);
    //char str[2];
    double mass = 0;
    for (int i = 0; i < strl; i++) {
        //*(str+0) = p[i];
        //*(str+1) = '\0';
        mass = mass + aa_table[p[i] - 42]; //AminoAcidMass(p[i]); //getAA(str);
        //printf("Mass is %f.\n", mass);
    }
    return mass + H2O;
}

