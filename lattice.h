#ifndef _LATTICE_H
#define _LATTICE_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"

/* -------------------------------------------------------------------- */
/**
 * Inline functions
 */
#define put_to(mat, i, j, val) mpz_set(mat[i][j+1].c, val)
#define smult_lattice(mat, i, j, factor) mpz_mul(mat[i][j+1].c, mat[i][j+1].c, factor)
#define get_entry(mat, i, j) mat[i][j+1].c

extern void handle_signals(lattice_t *lattice, DOUBLE **R);
extern void stop_program_sig(int sig);

/* Basic subroutines */
extern void debug_print(char *m, int l);
extern void print_lattice(lattice_t *lattice);
extern void print_lattice_stat(lattice_t *lattice, DOUBLE **R);
extern void dump_lattice(lattice_t *lattice);
extern void load_lattice(lattice_t *lattice, char *fname);

extern long gcd(long n1, long n2);
extern void coeffinit(coeff_t *v, int z);

extern DOUBLE scalarproductlfp (coeff_t *v, coeff_t *w);
extern DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w , int n);

extern int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N,  DOUBLE ***bs, int s, int z);
extern int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s);
extern double orthogonality_defect(lattice_t *lattice, DOUBLE **R, int s, int z);
extern double log_potential(DOUBLE **R, int s, int z);

extern int log2mpz(mpz_t number);
extern int get_bit_size(lattice_t *lattice);

extern void shufflelattice(lattice_t *lattice);
#endif
