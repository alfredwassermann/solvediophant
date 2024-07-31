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

extern DOUBLE exhaustive_enumeration(lattice_t *lattice);

extern DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max);
extern DOUBLE compute_w2(DOUBLE *w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
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
