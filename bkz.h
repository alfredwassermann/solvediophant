#ifndef _BKZ_H
#define _BKZ_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"

#define ENUM_BLOCK  0
#define ENUM_LDS_FULL  1
#define ENUM_LDS_FULL2 2
#define ENUM_LDS_BLOCK 3

/* -------------------------------------------------------------------- */
/**
 * Helper arrays for the function enumerate in bkz
 */
typedef struct {
    DOUBLE *c;
    DOUBLE *y;
    long *delta;
    long *d;
    long *v;
    DOUBLE *u_loc;
} bkz_enum_t;

extern void allocate_bkz_enum(bkz_enum_t *bkz_enum, int s);
extern void free_bkz_enum(bkz_enum_t *bkz_enum);

extern DOUBLE bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, DOUBLE p,
    int enum_type, int max_tours,
    int (*solutiontest)(lattice_t *lattice, int k),
    int (*solutiontest_long)(lattice_t *lattice, int k));

extern DOUBLE enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block,
                DOUBLE improve_by, DOUBLE p, bkz_enum_t *bkz_enum);
extern DOUBLE lds_enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block,
                DOUBLE improve_by, DOUBLE p, bkz_enum_t *bkz_enum);

extern void insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv);
extern void insert_vector_long(lattice_t *lattice, long *u, int start, int end, int z);

extern DOUBLE GH(DOUBLE **R, int low, int up);
extern void hoerner(DOUBLE **R, int low, int up, double p, DOUBLE *eta);
extern DOUBLE set_prune_const(DOUBLE **R, int low, int up, int prune_type, DOUBLE p);

#endif
