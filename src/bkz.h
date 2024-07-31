/**
Copyright 2024 Alfred Wassermann

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _BKZ_H
#define _BKZ_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"

#define ENUM_BLOCK     101
#define ENUM_LDS_FULL  201
#define ENUM_LDS_FULL2 202
#define ENUM_LDS_BLOCK 203

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
