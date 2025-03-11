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

extern double exhaustive_enumeration(lattice_t *lattice);

extern double compute_y(double **mu_trans, double *us, int level, int level_max);
extern double compute_w2(double *w, double **bd, double alpha, int level, int rows);
extern double compute_w(double *w, double *w1, double **bd, double alpha, int level, int rows);
extern double compute_w2float(float *w, float **bd, double alpha, int level, int rows);
extern double compute_wfloat(float *w, float *w1, float **bd, double alpha, int level, int rows);
extern void gramschmidt(lattice_t *lattice, int columns, int rows, double **mu, double **bd, double *c);
extern void givens(lattice_t *lattice, int columns, int rows, double **mu, double **bd, double *c);
extern void inverse(double **mu, double **muinv, int columns);
extern int final_test(double *v, int rows, double Fq, double *us, lattice_t *lattice);
// extern int prune(double *w, double cs, int rows, double Fqeps);
extern int prune_only_zeros(lattice_t *lattice, double *w, double *w1,
                int level, int rows, double Fq,
                double **bd, double y, int columns);
extern int prune_only_zeros_float(lattice_t *lattice, float *w, float *w1,
    int level, int rows, double Fq,
    float **bd, double y, int columns);
    
extern int print_solution(lattice_t *lattice, double *w, int rows, double Fq, double *us, int columns);
extern void print_num_solutions(long num_solutions);

#endif
