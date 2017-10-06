#ifndef _DATASTRUCT_H
#define _DATASTRUCT_H
#include <gmp.h>

#define SQRT sqrt
#define DOUBLE double

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
    int iterate;
    int iterate_no;

    int silent;
    long stop_after_solutions;
    long stop_after_loops;
    int print_ntl;

    struct _bkz {
        int beta;
        DOUBLE p;
    } bkz;

    struct _exhaustive_enum {
        int lds;
        int lds_k_max;
    } exhaustive_enum;
    
} lll_params_t;

/**
 * @name decomp_t
 * Helper arrays for Gram-Schmidt dcomposition and others
 */
typedef struct {
    DOUBLE *bd_1norm;
    int *first_nonzero, *first_nonzero_in_column, *firstp;
    int bit_size;

    DOUBLE *N;
    DOUBLE *c;
    DOUBLE **bd; 
    DOUBLE **mu, **mu_trans;

    DOUBLE Fd, Fq, Fqeps;
    DOUBLE tmp;
    coeff_t *swap_vec;

    DOUBLE stepWidth; // = 0.0;
    DOUBLE old_coeff;

    DOUBLE *fipo;
    DOUBLE **dual_basis;
    DOUBLE *dual_bound;
} decomp_t;

/**
 * @name lattice_t
 */
typedef struct {
    int num_rows;
    int num_cols;
    coeff_t **basis;
    long **basis_long;
    coeff_t *swap;
    long *swap_long;

    mpz_t matrix_factor;
    mpz_t max_norm;
    int cut_after;
    int free_RHS;

    int work_on_long;

    decomp_t *decomp;
    
    lll_params_t LLL_params;
    
} lattice_t;

/**
 * Needed for swaping in bkz, see also dio2.c
 */
#define ADDITIONAL_COLS 1

#endif
