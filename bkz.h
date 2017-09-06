#ifndef _BKZ_H
#define _BKZ_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"

/* -------------------------------------------------------------------- */

extern DOUBLE bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, DOUBLE p,
    int (*solutiontest)(lattice_t *lattice, int k),
    int (*solutiontest_long)(lattice_t *lattice, int k));

extern DOUBLE enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block, DOUBLE improve_by, DOUBLE p);
extern DOUBLE lds_enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block, DOUBLE improve_by, DOUBLE p);

extern void insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv);
extern void insert_vector_long(lattice_t *lattice, long *u, int start, int end, int z);

extern DOUBLE GH(DOUBLE **R, int low, int up);
extern void hoerner(DOUBLE **R, int low, int up, double p, DOUBLE *eta);
extern DOUBLE set_prune_const(DOUBLE **R, int low, int up, int prune_type, DOUBLE p);

#endif
