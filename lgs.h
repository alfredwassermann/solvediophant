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

    int num_original_cols;
    int *original_cols;
} lgs_t;

#define ZLENGTH 16000

extern void lgs_allocate_mem(lgs_t *LGS);
extern void lgs_free_mem(lgs_t *LGS);
extern void read_upper_bounds(char *file_name, lgs_t *LGS);
extern void read_selected_cols(char *file_name, lgs_t *LGS);
extern void read_linear_system(FILE *txt, lgs_t *LGS);
extern void incorrect_input_file();

#endif
