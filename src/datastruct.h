#ifndef _DATASTRUCT_H
#define _DATASTRUCT_H
#include <stdbool.h>
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
    int type;
    int iterate_no;

    int silent;
    long stop_after_solutions;
    long stop_after_loops;
    int print_ntl;

    struct _lll {
        DOUBLE delta_low;
        DOUBLE delta_med;
        DOUBLE delta_high;
        DOUBLE delta_higher;
    } lll;

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
    int bit_size;

    DOUBLE *N;
    DOUBLE *c;
    DOUBLE **bd;
    DOUBLE **mu, **mu_trans;

    DOUBLE *h_beta;
    DOUBLE **H;
    DOUBLE **R;

    DOUBLE Fd, Fq, Fqeps;
    DOUBLE tmp;
    mpz_t *swap_vec;

    DOUBLE stepWidth; // = 0.0;
    DOUBLE old_coeff;

    DOUBLE *fipo;
    DOUBLE **dual_basis;
    DOUBLE *dual_bound;

    DOUBLE *bd_1norm;

    int *first_nonzero;
    int *first_nonzero_in_column;
    int *firstp;
} decomp_t;

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
    long **basis_long;
    mpz_t *swap;
    long *swap_long;

    mpz_t matrix_factor;
    mpz_t max_norm;
    mpz_t max_norm_initial;

    mpz_t *upperbounds;
    mpz_t upperbounds_max;
    //mpz_t upfac;
    mpz_t max_up;

    int *original_cols;
    int no_original_cols;

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
