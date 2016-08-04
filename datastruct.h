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
        int p;
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
    //coeff_t **basis_s;
    //
    long **basis_long;
    coeff_t *swap;

    mpz_t matrix_factor;
    mpz_t max_norm;
    int cut_after;
    int free_RHS;

    lll_params_t LLL_params;
} lattice_t;

#define ADDITIONAL_COLS 20

#endif
