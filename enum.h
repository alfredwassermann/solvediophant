#ifndef _ENUM_H
#define _ENUM_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"

/* -------------------------------------------------------------------- */
typedef struct {
    mpz_t u;
    mpz_t s;
    mpz_t upfac;
    FILE* fp;
} solution_t;

extern DOUBLE explicit_enumeration(lattice_t *lattice);

extern DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max);
extern void compute_w2(DOUBLE *w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern DOUBLE compute_w(DOUBLE *w, DOUBLE *w1, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void gramschmidt(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c);
extern void givens(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c);
extern void inverse(DOUBLE **mu, DOUBLE **muinv, int columns);
extern int final_test(DOUBLE *v, int rows, DOUBLE Fq, DOUBLE *us, lattice_t *lattice);
extern int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps);
extern int prune_only_zeros(lattice_t *lattice, DOUBLE *w, DOUBLE *w1,
                int level, int rows, DOUBLE Fq,
                DOUBLE **bd, DOUBLE y, int columns);

extern int print_solution(lattice_t *lattice, DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns);
extern void print_num_solutions(long num_solutions);

#endif
