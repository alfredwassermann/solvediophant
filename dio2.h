#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h>
#include "lgs.h"

#define SQRT sqrt
#define DOUBLE double

/**
 * Definition of the lattice data structures
*/
typedef struct {
    mpz_t c;
    int p;
} coeff_t;

typedef struct {
    mpz_t scalelastlinefactor;
    int iterate;
    int iterate_no;
    struct Bkz {
        int beta;
        int p;
    } bkz;

    int silent;
    long stop_after_solutions;
    long stop_after_loops;
    int print_ntl;
} lll_params_t;

typedef struct {
    int num_rows;
    int num_cols;
    coeff_t **basis;
    //coeff_t **basis_s;
    coeff_t *swap;

    mpz_t matrix_factor;
    mpz_t max_norm;
    int cut_after;
    int free_RHS;

    lll_params_t LLL_params;
} lattice_t;
#define ADDITIONAL_COLS 20

/* Global variable used in stop_program */
int SILENT;
int PRINT_REQUIRED;
int DUMP_REQUIRED;

/* -------------------------------------------------------------------- */

extern long diophant(lgs_t *LGS, lattice_t *lattice, FILE* solfile, int restart, char *restart_filename);

extern long nosolutions;

extern void handle_signals(lattice_t *lattice);
extern void stop_program_sig(int sig);
extern void print_lattice_sig(int sig);
extern void dump_lattice_sig(int sig);

/* Basic subroutines */
extern void print_num_solutions(long num_solutions);
extern void debug_print(char *m, int l);
extern void print_lattice(lattice_t *lattice);
extern void dump_lattice(lattice_t *lattice);
extern void load_lattice(lattice_t *lattice, char *fname);

extern long gcd(long n1, long n2);
extern void coeffinit(coeff_t *v, int z);
extern int cutlattice(lattice_t *lattice);
extern int solutiontest(lattice_t *lattice, int position);

extern DOUBLE scalarproductlfp (coeff_t *v, coeff_t *w);
extern DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w , int n);

extern void check_precision(coeff_t *b, DOUBLE *R, int z, int k);

extern int householder_column(coeff_t **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z, int bit_size);
extern void size_reduction(coeff_t **b, DOUBLE **mu, mpz_t musvl, DOUBLE mus, int k, int j);

extern int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N,  DOUBLE ***bs, int s, int z);
extern int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s);
extern double orthogonality_defect(lattice_t *lattice, DOUBLE **R, int s, int z);
extern double log_potential(DOUBLE **R, int s, int z);

extern int log2mpz(mpz_t number);
extern int get_bit_size(lattice_t *lattice);

extern void lll(lattice_t *lattice, int s, int z, DOUBLE quality, int deepinsert_blocksize);
extern DOUBLE iteratedlll(lattice_t *lattice, int s, int z, int no_iterates, DOUBLE quality, int deepinsert_blocksize);
extern DOUBLE block_reduce(lattice_t *lattice, int s, int z, int block_size, DOUBLE quality, int deepinsert_blocksize);

extern DOUBLE bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, int p);
extern DOUBLE enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block, DOUBLE improve_by, int p);
extern DOUBLE sample(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block);
extern DOUBLE explicit_enumeration(lattice_t *lattice, int columns, int rows);

extern DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max);
extern void compute_w2(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void gramschmidt(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq);
extern void givens(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq);
extern void inverse(DOUBLE **mu, DOUBLE **muinv, int columns);
extern int final_test(DOUBLE *v, int rows, DOUBLE Fq, DOUBLE *us, lattice_t *lattice, int bit_size);
extern int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps);
extern int prune_only_zeros(DOUBLE **w, int level, int rows, DOUBLE Fq,
                     int *first_nonzero_in_column, int *firstp,
                     DOUBLE **bd, DOUBLE *y, DOUBLE *us, int columns);
extern int print_solution(lattice_t *lattice, DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns);
extern void shufflelattice(lattice_t *lattice);
extern DOUBLE GH(DOUBLE **R, int low, int up);
extern void hoerner(DOUBLE **R, int low, int up, double p, DOUBLE *eta);
extern DOUBLE set_prune_const(DOUBLE **R, int low, int up, int prune_type, DOUBLE p);
extern void print_NTL_lattice(lattice_t *lattice);
#endif
