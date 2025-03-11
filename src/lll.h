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

#ifndef _LLL_H
#define _LLL_H
#include <gmp.h>
#include "const.h"
#include "datastruct.h"
#include "lgs.h"

extern int lllH(lattice_t *lattice, double **R, double *beta, double **H,
            int start, int low, int up, int z,
            double delta, int reduction_type,
            int bit_size,
            int word_len,
            int (*solutiontest)(lattice_t *lattice, int k));

extern double householder_column(mpz_t **b, double **R, double **H, double *beta, int k, int z, int bit_size);
extern double householder_column_long(long **b, double **R, double **H, double *beta, int k, int z, int bit_size);

extern void size_reduction(mpz_t **b, double **R, mpz_t musvl, double mus, int k, int j, int z);
extern void size_reduction_long(long **b, double **R, long musl, double mus, int k, int j, int z);

extern void check_precision(mpz_t *b, double *R, int z, int k);

#endif
