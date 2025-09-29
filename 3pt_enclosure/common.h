#ifndef COMMON_H_
#define COMMON_H_


/* Program     : loop                                                       */
/* Compilation : make                                                       */
/* Structure constants in O(n) model: dense/dilute phase on square lattice. */
/* Periodic boundary conditions. Weights as in Warnaar-Nienhuis-Seaton.     */
/* Multiplied by -2 sin(2l)cos(3l). Take u=-3l/2+3Pi/4 (isotropic point).   */
/* Dense regime: 0 < l <= Pi/4. Dilute regime: Pi/4 <= l < Pi/2.            */
/* Version 01/10/2013. */
/* Allow bottom strings to loop around middle point. Version 03/04/2024.    */
/* Same but with diagonal operator in the middle. Version 16/04/2024.       */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Constants to be altered before each run ******************************** */

#define n                  4 /* Strip width/number of strands (any parity) */
#define M              (20*n) /* Half length */
#define hash_size    8000001

#define kk1                 2 /* Number of leg insertions (any parity)      */
#define kk3                 2
#define n2                2.0 /* Constant loop weight for diagonal operator */

// #define RESUME           0.250 /* Weight lambda_over_pi where to resume uncomplete computation (comment out to compute all).  */

// #define DEBUG
// #define VERBOSE

/* Other constants ******************************************************** */

/* For the usual version without looping strings */
#define EMPTY    0
#define BOTTOM   1
#define MIDDLE   2
#define ORDINARY (kk1+1)

#define PI  3.1415926535897932385
#define N            (n+2)   /* Maximal number of dangling ends */
#define hash_modulo    251
#define size             N

/* Global variables ******************************************************* */

double n_loop;
double rho[10]; /* Boltzmann weights */
double exponent;
int min_strings;
int target[kk1+1];

typedef struct entry         /* One entry in the hash table */
{
  double weight;
  struct entry *p_next;
  char key[size];
} T_entry;

T_entry* p_hash[2][hash_size]; /* The hash tables: We need 5 of them */
double weight;
char key[size];
int data_position;

/* These are the function declarations ************************************ */

double compute_Z(int k1,int k2,int k3);
double compute_Z_special(int k1,int k2,int k3);
double compute_Z_seam(int lower_seam,int upper_seam);
void transfer_row();
void transfer_row_special(int seam);
void insert_aux_space();
void insert_aux_space_special(int seam);
void R_matrix(int pos);
void R_matrix_special(int pos);
void remove_aux_space();
void remove_aux_space_special();
double getloopweight(unsigned char k);
double vector_norm(int data_pos);
void normalise_vector(int data_pos);
void insert(int data_pos); /* Insert in hash table */
void insert_special(int data_pos); /* Insert in hash table */
double get_weight(int data_pos);
void write_hash(int data_pos);
void write_hash_statistics(int data_pos);
void write(char *number);
void Strcpy(char *key1,char *key2);
int Strcmp(char *key1,char *key2);
void eliminate(T_entry* p_entry);
void empty_hash(int data_pos);
void init_hash(void);

#endif // end COMMON_H_
