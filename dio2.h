#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h>
#include "gls.h"

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
    int num_rows;
    int num_cols;
    coeff_t **basis;
    coeff_t *swap;

    mpz_t matrix_factor;
    mpz_t max_norm;
    int cut_after;
    int free_RHS;    
} lattice_t;

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

} lll_params_t;

/* -------------------------------------------------------------------- */

extern long diophant(gls_t *GLS, lll_params_t *LLL_params, FILE* solfile);

extern long nosolutions;

extern void stop_program(int sig);
extern void show_lattice(int sig);

/* Basic subroutines */
extern void print_num_solutions(long num_solutions);
extern void debug_print(char *m, int l);
extern void print_lattice();
extern long gcd(long n1, long n2);
extern void coeffinit(coeff_t *v, int z);
extern int cutlattice();
extern int solutiontest(int position);

extern DOUBLE scalarproductlfp (coeff_t *v, coeff_t *w);
extern DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w , int n);

extern void check_precision(coeff_t *b, DOUBLE *R, int z, int k);

extern int householder_column(coeff_t **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z);
extern void size_reduction(coeff_t **b, DOUBLE **mu, mpz_t musvl, DOUBLE mus, int k, int j);

extern int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N,  DOUBLE ***bs, int s, int z);
extern int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s);
extern double orthogonality_defect(coeff_t **lattice, DOUBLE **R, int s, int z);
extern double log_potential(DOUBLE **R, int s, int z);

extern void lll(coeff_t **b, int s, int z, DOUBLE quality, int deepinsert_blocksize);
extern DOUBLE iteratedlll(coeff_t **b, int s, int z, int no_iterates, DOUBLE quality, int deepinsert_blocksize);
extern DOUBLE bkz(coeff_t **b, int s, int z, DOUBLE delta, int beta, int p);
extern DOUBLE enumerate(DOUBLE **mu, DOUBLE *c, long *u, int s, int start_block, int end_block, int p);
extern DOUBLE explicit_enumeration(coeff_t **lattice, int columns, int rows);

extern DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max);
extern void compute_w2(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void gramschmidt(coeff_t **lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq);
extern void givens(coeff_t **lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq);
extern void inverse(DOUBLE **mu, DOUBLE **muinv, int columns);
extern int exacttest(DOUBLE *v, int rows, DOUBLE Fq);
extern int prune0(DOUBLE li, DOUBLE re);
extern int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps);
extern int prune_only_zeros(DOUBLE **w, int level, int rows, DOUBLE Fq,
                     int *first_nonzero_in_column, int *firstp,
                     DOUBLE **bd, DOUBLE *y, DOUBLE *us, int columns);
extern int print_solution(DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns);
extern void shufflelattice();
#endif
