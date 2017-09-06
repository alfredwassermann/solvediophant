#ifndef _DUALBKZ_H
#define _DUALBKZ_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"

extern DOUBLE dual_bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, DOUBLE p,
    int (*solutiontest)(lattice_t *lattice, int k));
extern DOUBLE self_dual_bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, DOUBLE p,
    int (*solutiontest)(lattice_t *lattice, int k));
extern void dual_insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv);
