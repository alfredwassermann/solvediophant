#ifndef _DATASTRUCT_H
#define _DATASTRUCT_H
#include <gmp.h>

#define SQRT sqrt
#define DOUBLE double

/**
 * Definition of the lattice data structures
*/
typedef struct {
    mpz_t c;
    int p;
} coeff_t;

typedef struct {
    mpz_t scalelastlinefactor;
    int iterate;
    int iterate_no;
    struct Bkz {
        int beta;
        DOUBLE p;
    } bkz;

    int silent;
    long stop_after_solutions;
    long stop_after_loops;
    int print_ntl;
} lll_params_t;

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

    lll_params_t LLL_params;
} lattice_t;

#define ADDITIONAL_COLS 1

#endif
