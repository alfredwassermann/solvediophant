#include <stdlib.h>
#include "gls.h"

/**
 * Allocate memnory for input matrix
 */
void gls_allocate_mem(gls_t *GLS) {
    int i, j;

    GLS->matrix = (mpz_t**)calloc(GLS->num_rows, sizeof(mpz_t*));

    for (j = 0; j < GLS->num_rows; j++) {
        GLS->matrix[j] = (mpz_t*)calloc(GLS->num_cols, sizeof(mpz_t));
        for (i = 0; i < GLS->num_cols; i++)
           mpz_init(GLS->matrix[j][i]);
    }

    GLS->rhs = (mpz_t*)calloc(GLS->num_rows, sizeof(mpz_t));
    for (i = 0; i < GLS->num_rows; i++)
       mpz_init(GLS->rhs[i]);
}

void gls_free_mem(gls_t *GLS) {
    int i, j;

    for (j = 0; j < GLS->num_rows; j++) {
        for (i = 0; i < GLS->num_cols; i++)
			mpz_clear(GLS->matrix[j][i]);
        free(GLS->matrix[j]);
    }
    free(GLS->matrix);

    for (i = 0; i < GLS->num_rows; i++)
		mpz_clear(GLS->rhs[i]);
    free(GLS->rhs);

    if (GLS->upperbounds != NULL) {
        for (i = 0; i < GLS->num_cols; i++) {
            mpz_clear(GLS->upperbounds[i]);
        }
        free(GLS->upperbounds);
    }
}
