#ifndef _GLS_H
#define _GLS_H
#include <stdio.h>
#include <gmp.h>

typedef struct {
    int num_rows;
    int num_cols;

    mpz_t **matrix;
    mpz_t *rhs;
    mpz_t *upperbounds;

    int num_boundedvars;

    int num_original_cols;
    int *original_cols;
} gls_t;

#define ZLENGTH 16000

extern void gls_allocate_mem(gls_t *GLS);
extern void gls_free_mem(gls_t *GLS);
extern void read_upper_bounds(char *file_name, gls_t *GLS);
extern void read_selected_cols(char *file_name, gls_t *GLS);
extern void read_linear_system(FILE *txt, gls_t *GLS);
extern void incorrect_input_file();

#endif
