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

#ifndef _DATASTRUCT_H
#define _DATASTRUCT_H
#include <stdbool.h>
#include <gmp.h>

#define SQRT sqrt
#define double double

/**
 * @name coeff_t
 * Optimize sparse integer vector operations
 */
typedef struct {
    mpz_t c;
    int p;
} coeff_t;

/**
 *  @name lll_params_t
 */
typedef struct {
    mpz_t scalelastlinefactor;
    int type;
    int iterate_no;

    int silent;
    long stop_after_solutions;
    long stop_after_loops;
    long report_interval;
    int print_ntl;  
    int kernel;


    struct _lll {
        double delta_low;
        double delta_med;
        double delta_high;
        double delta_higher;
    } lll;

    struct _bkz {
        int beta;
        double p;
    } bkz;

    struct _exhaustive_enum {
        int lds;
        int lds_k_max;
    } exhaustive_enum;

} lll_params_t;

/**
 * @name decomp_t
 * Helper arrays for Gram-Schmidt decomposition and others
 */
typedef struct {
    int bit_size;

    double *N;
    double *c;
    double **bd;
    double **mu, **mu_trans;
    double *bdMemory, *muMemory; // Just for 32byte aligning

    double *h_beta;
    double **H;
    double **R;

    double Fd, Fq, Fqeps;
    double tmp;
    mpz_t *swap_vec;

    double stepWidth; // = 0.0;
    double old_coeff;

    double *fipo;
    double **dual_basis;
    double *dual_bound;

    double *bd_1norm;

    int *first_nonzero;
    int *first_nonzero_in_column;
    int *firstp;

    float **bdfloat;
    float **dummy1, **dummy2, **dummy3;  // Just for 32byte aligning

} decomp_t; //  __attribute__ ((aligned (32))) ;

/**
 * @name lattice_t
 */
typedef struct {
    int num_rows;
    int num_cols;
    int lgs_rows;
    int lgs_cols;
    int lgs_rank;

    int num_boundedvars;
    bool is_zero_one; // true, if 0/1 variables
    bool work_on_long;
    int free_RHS;
    int cut_after;

    mpz_t **basis;
    mpz_t *swap;
    long **basis_long;
    long *swap_long;

    mpz_t matrix_factor;
    mpz_t max_norm;
    mpz_t max_norm_initial;

    mpz_t *upperbounds;
    mpz_t upperbounds_max;
    //mpz_t upfac;
    mpz_t max_up;

    int *original_cols;
    int num_cols_org;
    int num_rows_org;

    long nom;
    long denom;

    decomp_t decomp;

    lll_params_t LLL_params;

} lattice_t;

typedef struct {
    long num_solutions;
} stat_t;

/**
 * Needed for swapping in bkz, see also dio2.c
 */
#define ADDITIONAL_COLS 1

#endif
