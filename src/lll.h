#ifndef _LLL_H
#define _LLL_H
#include <gmp.h>
#include "const.h"
#include "datastruct.h"
#include "lgs.h"

extern int lllH(lattice_t *lattice, DOUBLE **R, DOUBLE *beta, DOUBLE **H,
            int start, int low, int up, int z,
            DOUBLE delta, int reduction_type,
            int bit_size,
            int (*solutiontest)(lattice_t *lattice, int k));

extern int lllH_long(lattice_t *lattice, DOUBLE **R, DOUBLE *beta, DOUBLE **H,
            int start, int low, int up, int z,
            DOUBLE delta, int reduction_type,
            int bit_size,
            int (*solutiontest)(lattice_t *lattice, int k));

extern int householder_column(coeff_t **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z, int bit_size);
extern int householder_column_long(long **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z, int bit_size);

extern void size_reduction(coeff_t **b, DOUBLE **mu, mpz_t musvl, DOUBLE mus, int k, int j);
extern void size_reduction_long(long **b, DOUBLE **mu, long musl, DOUBLE mus, int k, int j, int z);

extern void check_precision(coeff_t *b, DOUBLE *R, int z, int k);

#endif
