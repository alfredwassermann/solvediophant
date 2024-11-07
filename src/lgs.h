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

#ifndef _LGS_H
#define _LGS_H
#include <stdio.h>
#include <gmp.h>

typedef struct {
    int num_rows;
    int num_cols;

    mpz_t **matrix;
    mpz_t *rhs;
    mpz_t *upperbounds;

    int num_boundedvars;
    int rank;

    int num_original_cols;
    int *original_cols;
} lgs_t;

#define ZLENGTH 16000

extern void alloc_lgs(lgs_t *LGS);
extern void free_lgs(lgs_t *LGS);
extern void read_upper_bounds(FILE *txt, char *zeile, lgs_t *LGS);
extern void read_selected_cols(FILE *txt, lgs_t *LGS);
extern void read_linear_system(FILE *txt, lgs_t *LGS);
extern int  rank(lgs_t *LGS, long p);
extern int  preprocess(lgs_t *LGS);
extern void incorrect_input_file();

#endif
