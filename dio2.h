#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"

/* -------------------------------------------------------------------- */

extern long diophant(lgs_t *LGS, lattice_t *lattice, FILE* solfile, int restart, char *restart_filename);

extern long nosolutions;

/* Basic subroutines */
extern int cutlattice(lattice_t *lattice);
extern int solutiontest(lattice_t *lattice, int position);
extern int solutiontest_long(lattice_t *lattice, int position);

extern void lll(lattice_t *lattice, int s, int z, DOUBLE quality, int reduction_type);
extern DOUBLE iteratedlll(lattice_t *lattice, int s, int z, int no_iterates, DOUBLE quality, int reduction_type);
extern DOUBLE block_reduce(lattice_t *lattice, int s, int z, int block_size, DOUBLE quality, int reduction_type);

extern DOUBLE bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, int p);
extern DOUBLE dual_bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, int p);
extern DOUBLE enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block, DOUBLE improve_by, int p);
extern DOUBLE dual_enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block, DOUBLE improve_by, int p);
extern void insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv);
extern void insert_vector_long(lattice_t *lattice, long *u, int start, int end, int z);
extern void dual_insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv);

extern DOUBLE sample(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block);
extern DOUBLE explicit_enumeration(lattice_t *lattice, int columns, int rows);

extern DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max);
extern void compute_w2(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void gramschmidt(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c);
extern void givens(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c);
extern void inverse(DOUBLE **mu, DOUBLE **muinv, int columns);
extern int final_test(DOUBLE *v, int rows, DOUBLE Fq, DOUBLE *us, lattice_t *lattice, int bit_size);
extern int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps);
extern int prune_only_zeros(DOUBLE **w, int level, int rows, DOUBLE Fq,
                     int *first_nonzero_in_column, int *firstp,
                     DOUBLE **bd, DOUBLE *y, DOUBLE *us, int columns);

extern int print_solution(lattice_t *lattice, DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns);

extern DOUBLE GH(DOUBLE **R, int low, int up);
extern void hoerner(DOUBLE **R, int low, int up, double p, DOUBLE *eta);
extern DOUBLE set_prune_const(DOUBLE **R, int low, int up, int prune_type, DOUBLE p);
extern void print_NTL_lattice(lattice_t *lattice);
extern void print_num_solutions(long num_solutions);

extern void print_lattice_sig(int sig);
extern void dump_lattice_sig(int sig);
#endif
