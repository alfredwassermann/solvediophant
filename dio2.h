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

extern DOUBLE explicit_enumeration(lattice_t *lattice, int columns, int rows);

extern DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max);
extern void compute_w2(DOUBLE *w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern DOUBLE compute_w(DOUBLE *w, DOUBLE *w1, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void gramschmidt(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c);
extern void givens(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c);
extern void inverse(DOUBLE **mu, DOUBLE **muinv, int columns);
extern int final_test(DOUBLE *v, int rows, DOUBLE Fq, DOUBLE *us, lattice_t *lattice, int bit_size);
extern int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps);
extern int prune_only_zeros(lattice_t *lattice, DOUBLE *w, DOUBLE *w1,
                int level, int rows, DOUBLE Fq,
                int *first_nonzero_in_column, int *firstp,
                DOUBLE **bd, DOUBLE y, int columns);

extern int print_solution(lattice_t *lattice, DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns);

extern void print_NTL_lattice(lattice_t *lattice);
extern void print_num_solutions(long num_solutions);

extern void print_lattice_sig(int sig);
extern void dump_lattice_sig(int sig);
#endif
