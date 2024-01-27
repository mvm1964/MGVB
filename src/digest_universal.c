/*
Digester program for the LC-MS/MS package    
Reads from fasta file and produces a txt file
with tab-delimited sequence, mass and header  
The file is then used by toSQL to    
import into sqlite3 database, which is used  
by the scorrer program to compute scores.

Also produces reversed peptides for decoy searches but this might be
outsourced later to the scorer program and don on the fly to save space
and optimise search.
    
Usage: digest <fasta file> <peptide file> <n missed>
                                             
Metodi V. Metodiev, University of Essex, 2019 

Recursive algorithm for digesting proteins with n missed cleavages

MVFTRGPPLTTEWSKTTYEWSQRTTTGGVCA

1. Remove N-terminal M if requested
2. Scan string until first cleavage site and cut. Save peptide to file if conditions met.
3. Scan strings until second site found. Save string from start to this site as a peptide with 1 missed site.
4. Scan until next cleavage site. Save as above as peptide with 2 missed sites.
5. Repeat until a peptide with n missed sites is saved to file.
6. Now recursively execute 2-5 on a string starting just after the first cleavage site.

Update 29.06.2019
Start rethinking the database structure. Will test performance with more normalized database:
Create separate tables for:
   1. Proteins: should have a key, gene, sequence, and fasta header columns. The actual algorithm will be like this:
        - for each entry in the fasta file assign a number to be used as key. 
        - Write to a proteins file as the fasta file is being processedi.
        - Later bulk-load into sqlite3. 
   2. Peptides: should have a key, sequence and protein id (as foreign key) columns. Protein id is the key from Proteins 
      table. The entries should be non redundant so protein id column may contain multiple entries.
   3. Mod_peptides: should contain a key, peptide id (the key from Peptides table), mass, mod1, mod2, mod3...
      columns which are the corresponding names of the modifications containing numbers of such modifications. Should also 
      contain an entry for the unmodified peptide. This is the table that will be searched when there are modifications.
      This is likely to slow things though: instead of one sql query there will be multiples for each if the scorer remains as it is. 
      To speed it up the final report should be writen after all searches have finished so the scorer will first output sequence, 
      protein id, mass, score, ret time and scan number. Then the fasta header will be added for significant proteins at some
      point. It remains to be seen if keeping the sequence in a separate table will slow the search. If so it will have to be added
      to this table. It could be a small trade off in normalizing but an advantage in speed. The table should be indexed by peptide mass.
      How to handle decoy peptides: since decoy peptides will have the same modifications it would be more efficient to not compute multiset
      combinations for them but to just assign the same as for the true peptides. This means that the program that will generate the 
      mod_peptides need to be changed to only select non decoy peptides. One way to handle this is to generate decoy peptides on the fly
      in the scorer program which could work like that:
        - scorer receives matching entries from mod_peptides.
        - for each retrieves sequence from peptides and generates reverse sequence.
        - for true and reverse sequence calls a function which:  
          - Generates combinations of modified peptides using the mod numbers.
          - Computes score and outputs to file for each using the protein key and "decoy_" + protein key to link hits to proteins.

  Update 04.08.2020: change of strategy for decoy sequences. Reversing trypric peptides on the fly does not work. Too many high-scorring 
  decoy sequences. Will reverse protein sequences and digest them.

  Algorithm:
    1. Add a field "decoy" to the protein struct. Values will be 0 for not decoy and 1 for decoy
    2. In main funtion, reverse protein sequence, change decoy field to 1 and repeat the digesting   

  Update 16.08.2020: will modify it to accept a config file as argument and digest multiple fasta files:

    1. Use the code from scorer_15082020.c to parse the config file. Only fasta files are needed
    2. Process each file in a loop

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "pepMass.h"
#include <libgen.h>
#define maxLength 3500
#define lowerLimit 5
#define upperLimit 60

FILE *pepFl; // Might change this to not be a global pointer.

typedef struct Protein {
    char header[2048];
    char sequence[100000];
    char gene[2048]; 
    int decoy;
    int contam;
} protein;

// Empty peptide and protein.
protein emptyP = {};

// This function gets the gene name form fasta header from inside of protein struct.
void get_gene(protein *prot) {
    // Get Gene name.
    char gene[2048];
    strcpy(gene, prot->header);
    char *token;
    int i = 0;
    token = strstr(gene, "GN=");
    if (token != NULL) {
        while (!isspace(token[i])) i++;
        token[i] = '\0';
        strcpy(gene, token+3);
    }
    else strcpy(gene, "Unknown");
    strcpy(prot->gene, gene);
}


// This function reverses strings (copied it from the web). Will go to the new scorer.

char *strrev(char *str) {
      char *p1, *p2;

      if (! str || ! *str)
            return str;
      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}


// Recursive digest function. Will have to expand to include other enzymes.
/* Update 17.08.2020: other enzymes included but only C-terminal cutters. Will change to include N-term cutters too*/
void digest_with_missed(int prot_id, char *sequence, int n_missed, char* str1, int decoy, char *rule1, char *rule2, char site, int contam) {
    //char head[1024];
    //strcpy(head, header);     
    double pMass;
    int k = 0;
    int ind[n_missed + 1];
    ind[0] = strlen(sequence) - 1; // to prevent errors. Somehow on mac it works without it but on ubuntu it does not

    // Scan for cleavage sites and write to an index array. It considers C termini a cleavage site.
    for (int i = 0; i < strlen(sequence); i++) {
        if (strchr(rule1, sequence[i]) != NULL && strchr(rule2, sequence[i+1]) == NULL) { 
            ind[k] = i;
            k++;
            if (k > n_missed) break;    
        }
    }

    // Form peptide strings and write to file.
    for (int j = 0; j <= k -1; j++) {
        // Form sequence string. Note the way precision is used to determine length.

        if (site == 'C') sprintf(str1, "%.*s", ind[j] + 1, sequence);
        else sprintf(str1, "%.*s", ind[j], sequence);
        if (strlen(str1) <= upperLimit && strlen(str1) > lowerLimit) {
            pMass = pepMass(str1); // Calculate mass. Need to add Pr before searching. 
            if (fprintf(pepFl, "%s\t%lf\t%d\t%d\t%d\t%d\n", str1, pMass + Pr, prot_id, j, decoy, contam) == EOF) {
                printf("Error writing to file!\n");
                exit(1);
            }

            /*
            char revSeq[strlen(str1)];
            char s[strlen(str1)];
            strcpy(revSeq, str1);
            revSeq[strlen(str1)-1] = '\0';
            strcpy(s, strrev(revSeq));
            fprintf(pepFl, "%s\t%lf\t%s\t%s\t%d\n", strcat(s, str1 + strlen(str1) - 1), pMass + Pr, header, decoy, j);
            */
        }
    }
    if (sequence[ind[0] + 1] != '\0') {
        digest_with_missed(prot_id, sequence + ind[0] + 1, n_missed, str1, decoy, rule1, rule2, site, contam);
    }
    return;
} 

/* A middleware function. This is the first one to be called from main(). Since the digester function is recursive
   we need this one to stand in between, to extract needed data and pass it by reference to the digester..
*/
void digestMidware(protein *prot, int n_missed, char *str1, int prot_id, int decoy, char *rule1, char *rule2, char site, int contam) {
    // Removes N-terminal M.
    if (prot->sequence[0] == 'M') {
        prot->sequence[0] = '*'; 
        digest_with_missed(prot_id, prot->sequence, n_missed, str1, decoy, rule1, rule2, site, contam);
    }
    else digest_with_missed(prot_id, prot->sequence, n_missed, str1, decoy, rule1, rule2, site, contam); 
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Need to provide config file and peptide file names!\n");
        exit(1);
    }
    // Variable declarations.
    // this is causing valgrind errors. No, it is the missing not cut after value, now putting Z in config and no problems. 
    protein prot = emptyP;
    /*protein prot;
    strcpy(prot.header, "");
    strcpy(prot.sequence, "");
    strcpy(prot.gene, "");
    prot.decoy = 0;
    prot.contam = 0;*/

    char tmp_seq[100000];
    char fileStr[2048];
    char *dirName = dirname(argv[1]);

    // generate proteins file path and name
    sprintf(fileStr, "%s/proteins_from_fasta.txt", dirName);

    FILE *fptr;
    FILE *prtn;
    char str[maxLength];
    char str1[maxLength]; // To pass to digest function.

    // Open peptide file for writing.
    /*if ((fptr = fopen(argv[1], "r")) == NULL) {
        printf("Error opening the fasta file.\n");
        exit(1);
    }*/
    if ((pepFl = fopen(argv[2], "w")) == NULL) {
        printf("Error opening the peptide file.\n");
        exit(1);
    }
    
    int n_missed = 0;
    char digest_site = 'C';
    char digest_spec[16];
    char except_rule[16];
 
    printf("Preparing proteins. Please wait...\n");

    // read config file and parse fasta files names
    char *conf = argv[1];

    FILE *cnfg;
    char conf_str[1028];
    char fasta[5][2048];
    int fasta_cntr = 0;
    //int ptm_ctr = 0;

    // Open config file for reading
    if ((cnfg = fopen(conf, "r")) == NULL) {
        printf("Error opening config file!");
        exit(1);
    }
    // Read the config file line by line and parse it.
    while (fgets(conf_str, sizeof(conf_str), cnfg) != NULL) {

        if (strstr(conf_str, "_missed")) {
            sscanf(conf_str, "_missed\t%d", &n_missed); 
        }

        if (strstr(conf_str, "_digestSite")) {
            sscanf(conf_str, "_digestSite\t%c", &digest_site);
        }

        if (strstr(conf_str, "_cutAfter")) {
            sscanf(conf_str, "_cutAfter\t%s", digest_spec);
        }

        if (strstr(conf_str, "_notCutAfter")) {
            sscanf(conf_str, "_notCutAfter\t%s", except_rule);
        }

        if (strstr(conf_str, "_fasta")) {
            sscanf(conf_str, "_fasta\t%s", fasta[fasta_cntr]); // fasta should be crated in main with sufficient size
            fasta_cntr++;
        }
    }
    fclose(cnfg);
    
    // for debugging
    printf("Cutting %c-terminal of %s but not if followed by %s\n", digest_site, digest_spec, except_rule); 
    //exit(1);
    // Open a protein file to write to..    
    if ((prtn = fopen(fileStr, "w")) == NULL) {
        printf("Error opening the protein file.\n");
        exit(1);
    }
    int prot_id = 1; // this will be incremented and used as key for the proteins table.

    // Update 16.08.2020:Read fasta files in a loop line by line and parse protein sequences.
    for (int i = 0; i < fasta_cntr; i++) {
        printf("Processing %s. Please wait!\n", fasta[i]);
        if ((fptr = fopen(fasta[i], "r")) == NULL) {
            printf("Error opening the fasta file.\n");
            exit(1);
        }
        while (fgets(str, sizeof(str), fptr) != NULL) {
            str[strcspn(str, "\r\n")] = 0; // This is very important. Replaces returns and newlines with 0. 
            if (str[0] == '>') {
                // Digest previous prot struct if existing
                if (prot.sequence[0] != 0) {
                    // Write protein to file.
                    fprintf(prtn, "%d\t%s\t%s\t%s\t%d\t%d\n", prot_id, prot.gene, prot.header, prot.sequence, prot.decoy, prot.contam);
                    digestMidware(&prot, n_missed, str1, prot_id, prot.decoy, digest_spec, except_rule, digest_site, prot.contam);
                    
		    // for debugging
		    //printf("prot %d: %s\n", prot_id, prot.sequence);

		    prot_id++;
                }
                prot = emptyP;   
                strcpy(prot.header, str);
                get_gene(&prot);
                prot.decoy = 0;
                if (strstr(fasta[i], "contaminants")) prot.contam = 1;
                else prot.contam = 0;
            }
            else {
                strcat(prot.sequence, str);
            }
        }  
    
        // Add the last protein to the file and digest it
        fprintf(prtn, "%d\t%s\t%s\t%s\t%d\t%d\n", prot_id, prot.gene, prot.header, prot.sequence, prot.decoy, prot.contam);
        digestMidware(&prot, n_missed, str1, prot_id, prot.decoy,  digest_spec, except_rule, digest_site, prot.contam);
        prot_id++;
        prot = emptyP;
        fclose(fptr);
    }    

    // Now repeat everything with reverse sequences
    printf("Preparing decoy sequences...\n");

    // Read fasta file line by line and parse protein sequences.
    for (int i = 0; i < fasta_cntr; i++) {
        printf("Processing %s. Please wait!\n", fasta[i]);
        if ((fptr = fopen(fasta[i], "r")) == NULL) {
            printf("Error opening the fasta file.\n");
            exit(1);
        }
        while (fgets(str, sizeof(str), fptr) != NULL) {
            str[strcspn(str, "\r\n")] = 0; // This is very important. Replaces returns and newlines with 0.
            if (str[0] == '>') {
                // Digest previous prot struct if existing
                if (prot.sequence[0] != 0) {
                    // Write protein to file.
                    // Reverse sequence
                    strcpy(tmp_seq, strrev(prot.sequence));
                    strcpy(prot.sequence, tmp_seq);
                    fprintf(prtn, "%d\t%s\t%s\t%s\t%d\t%d\n", prot_id, prot.gene, prot.header, prot.sequence, prot.decoy, prot.contam);
                    digestMidware(&prot, n_missed, str1, prot_id, prot.decoy,  digest_spec, except_rule, digest_site, prot.contam);
                    prot_id++;
                }
                prot = emptyP;
                strcpy(prot.header, str);
                get_gene(&prot);
                prot.decoy = 1;
                if (strstr(fasta[i], "contaminants")) prot.contam = 1;
                else prot.contam = 0;
            }
            else {
                strcat(prot.sequence, str);
            }
        }
        // Add the last protein to the file and digest it
        strcpy(tmp_seq, strrev(prot.sequence));
        strcpy(prot.sequence, tmp_seq);
        fprintf(prtn, "%d\t%s\t%s\t%s\t%d\t%d\n", prot_id, prot.gene, prot.header, prot.sequence, prot.decoy, prot.contam);
        digestMidware(&prot, n_missed, str1, prot_id, prot.decoy,  digest_spec, except_rule, digest_site, prot.contam);
        //prot_id++;
        prot = emptyP;
        fclose(fptr);
    }

    //free(dirName);
    fclose(prtn);
    fclose(pepFl);
    printf("Digested %d proteins.\n", prot_id);       
    return 0;
} 
