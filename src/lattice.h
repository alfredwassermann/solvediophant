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

#ifndef _LATTICE_H
#define _LATTICE_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"

/* -------------------------------------------------------------------- */
/**
 * Inline functions
 */
#define put_to(mat, i, j, val) mpz_set(mat[i][j], val)
#define mult_by(mat, i, j, factor) mpz_mul(mat[i][j], mat[i][j], factor)
#define get_entry(mat, i, j) mat[i][j]

extern void handle_signals(lattice_t *lattice, DOUBLE **R);
extern void stop_program_sig(int sig);

/* Basic subroutines */
extern void lgs_to_lattice(lgs_t *LGS, lattice_t *lattice);

extern void debug_print(char *m, int l);
extern void print_lattice(lattice_t *lattice, FILE *stream);
extern void print_lattice_stat(lattice_t *lattice, DOUBLE **R);
extern void dump_lattice(lattice_t *lattice);
extern void load_lattice(lattice_t *lattice, char *fname);

extern long gcd(long n1, long n2);

extern DOUBLE dot_mpz(mpz_t *v, mpz_t *w, int z);

extern int alloc_decomp(lattice_t *lattice);
extern int free_decomp(decomp_t decomp); //DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s);
extern void free_lattice(lattice_t *lattice);
extern double orthogonality_defect(lattice_t *lattice, DOUBLE **R, int s, int z);
extern double log_potential(DOUBLE **R, int s, int z);

extern int log2mpz(mpz_t number);
extern int get_bit_size(lattice_t *lattice);

extern void shufflelattice(lattice_t *lattice);

extern void copy_lattice_to_long(lattice_t *lattice);
extern void copy_lattice_to_mpz(lattice_t *lattice);

extern void print_gsa(DOUBLE **R, int n, int start);

#endif
