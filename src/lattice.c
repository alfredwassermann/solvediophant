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

#include <signal.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <gmp.h>

#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"
#include "lll.h"
#include "arith.h"

#if defined(USE_BLAS)
    #define BLAS 1
    #include <cblas-openblas.h>
#elif defined(USE_BLAS_DEV)
    #define BLAS 1
    #include "common.h"
    #include "cblas.h"
#elif defined(USE_BLAS_OLD)
    #define BLAS 1
    #include <cblas.h>
#else
    #define BLAS 0
#endif


void debug_print(char *m, int l) {
    if (VERBOSE >= l) {
        printf("debug>> %s\n",m);
        fflush(stdout);
    }
    return;
}

/**
 * Print basis vectors as row vectors
 */
void print_lattice(lattice_t *lattice, FILE *stream) {
    int i, j;

    for (i = 0; i <= lattice->num_cols; i++) {
        for (j = 0; j < lattice->num_rows; j++) {
            if (lattice->work_on_long) {
                fprintf(stream, "%2ld ", lattice->basis_long[i][j]);
            } else {
                mpz_out_str(stream, 10, get_entry(lattice->basis, i, j));
                fprintf(stream, " ");
            }
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
    fflush(stream);
    return;
}

void dump_lattice(lattice_t *lattice) {
    int i, j;
    char fname[] = "dump_lattice.b";
    FILE* f = fopen(fname, "w");
    fprintf(f, "%d %d\n", lattice->num_rows, lattice->num_cols);
    fprintf(f, "%d\n", lattice->cut_after);
    fprintf(f, "%d\n", lattice->free_RHS);
    mpz_out_str(f, 10, lattice->matrix_factor); fprintf(f, "\n");
    mpz_out_str(f, 10, lattice->max_norm); fprintf(f, "\n");

    for (i = 0; i < lattice->num_cols; i++) {
        for (j = 0; j < lattice->num_rows; j++) {
            if (lattice->work_on_long) {
                fprintf(f, "%ld ", lattice->basis_long[i][j]);
            } else {
                mpz_out_str(f, 10, get_entry(lattice->basis, i, j));
                fprintf(f, " ");
            }
        }
        fprintf(f, "\n");
    }

    fclose(f);
    fprintf(stderr, "Lattice dumped to '%s'\n", fname);
    return;
}

void load_lattice(lattice_t *lattice, char *fname) {
    int i, j, res;

    fprintf(stderr, "LOAD lattice from file %s\n", fname);
    FILE* f = fopen(fname, "r");
    res = fscanf(f, "%d%d\n", &(lattice->num_rows), &(lattice->num_cols));
    if (res != 2) {
        exit(EXIT_ERR_INPUT);
    }
    res = fscanf(f, "%d\n", &(lattice->cut_after));
    if (res != 2) {
        exit(EXIT_ERR_INPUT);
    }
    res = fscanf(f, "%d\n", &(lattice->free_RHS));
        if (res != 1) {
        exit(EXIT_ERR_INPUT);
    }


    fprintf(stderr, "LOAD:  %d %d\n", lattice->num_rows, lattice->num_cols);

    res = mpz_inp_str(lattice->matrix_factor, f, 10);
    if (res == 0) { incorrect_input_file(); }
    res = mpz_inp_str(lattice->max_norm, f, 10);
    if (res == 0) { incorrect_input_file(); }

    for (i = 0; i < lattice->num_cols; i++) {
        for (j = 0; j < lattice->num_rows; j++) {
            res = mpz_inp_str(lattice->basis[i][j], f, 10);
            if (res == 0) {
                incorrect_input_file();
            }
        }
    }
    fclose(f);
    return;
}

/*
  Determine max bit size of lattice entries
 */
int log2mpz(mpz_t number) {
    int i = 1;
    mpz_t test, n;
    mpz_init_set_ui(test, 2);
    mpz_init(n);
    mpz_abs(n, number);

    while (mpz_cmp(test, n) < 0) {
        i++;
        mpz_mul_ui(test, test, 2);
    }
    mpz_clear(test);
    mpz_clear(n);
    return i;
}

int get_bit_size(lattice_t *lattice) {
    int i, j, log2_max = 0, log2_b;

    for (i = 0; i < lattice->num_cols; i++) {
        for (j = 0; j < lattice->num_rows; j++) {
            log2_b = log2mpz(get_entry(lattice->basis, i, j));
            if (log2_max < log2_b) {
                log2_max = log2_b;
            }
        }
    }

    return log2_max;
}

/**
 * LLL-subroutines
 */
DOUBLE dot_mpz (mpz_t *v, mpz_t *w, int z) {
    int i;
    mpz_t sum;

    mpz_init(sum);
    for (i = 0; i < z; ++i) {
        mpz_addmul(sum, v[i], w[i]);
    }

    return mpz_get_d(sum);
}

void alloc_basis (lattice_t *lattice) {
    int i, j;

    // Allocate memory for the gmp basis and the long basis
    lattice->basis = (mpz_t**)malloc((lattice->num_cols + ADDITIONAL_COLS) * sizeof(mpz_t*));
    lattice->basis_long = (long**)calloc(lattice->num_cols + ADDITIONAL_COLS, sizeof(long*));

    for (j = 0; j < lattice->num_cols + ADDITIONAL_COLS; j++) {
        lattice->basis[j] = (mpz_t*)malloc(((unsigned int)lattice->num_rows + 1) * sizeof(mpz_t));
        for (i = 0; i < lattice->num_rows + 1; i++) {
            mpz_init(lattice->basis[j][i]);
        }
        lattice->basis_long[j] = (long*)calloc((unsigned int)lattice->num_rows, sizeof(long));
    }

    // Allocate memory for swap vector
    lattice->swap = (mpz_t*)calloc((unsigned int)(lattice->num_rows + 1), sizeof(mpz_t));
    lattice->swap_long = (long*)calloc((unsigned int)(lattice->num_rows + 1), sizeof(long));
    for (i = 0; i < lattice->num_rows + 1; i++) {
        mpz_init(lattice->swap[i]);
    }
}

void free_lattice(lattice_t *lattice) {
    int i, j;

    // TODO: this is not exactly the allocated memory
    for (j = 0; j < lattice->num_cols + ADDITIONAL_COLS; j++) {
        for (i = 0; i < lattice->num_rows_org + 1; i++) { 
            mpz_clear(lattice->basis[j][i]);
        }
        free(lattice->basis[j]);
    }
    free(lattice->basis);

    for (j = 0; j < lattice->num_cols + ADDITIONAL_COLS; j++) {
        free(lattice->basis_long[j]);
    }
    free(lattice->basis_long);

    for (i = 0; i < lattice->num_rows_org + 1; i++) { mpz_clear(lattice->swap[i]); }
    free(lattice->swap);
    free(lattice->swap_long);

    mpz_clear(lattice->matrix_factor);
    mpz_clear(lattice->max_norm);
    mpz_clear(lattice->max_norm_initial);
    mpz_clear(lattice->max_up);
    mpz_clear(lattice->upperbounds_max);
    mpz_clear(lattice->LLL_params.scalelastlinefactor);

    if (lattice->upperbounds != NULL) {
        for (i = 0; i < lattice->lgs_cols /*lattice->num_cols_org*/; i++) {
            mpz_clear(lattice->upperbounds[i]);
        }
        free(lattice->upperbounds);
    }
    free(lattice->original_cols);

    free_decomp(lattice->decomp);
}

void handle_upperbounds(lgs_t *LGS, lattice_t *lattice) {
    int i;
    int lgs_cols = LGS->num_cols;

    // Handle upper bounds
    mpz_init_set_si(lattice->upperbounds_max, 1);

    lattice->is_zero_one = true;
    if (LGS->upperbounds == NULL) {
        fprintf(stderr, "No upper bounds: assume 0/1 variables \n"); fflush(stderr);
    } else {
        // Initialize the upper bounds with 1
        lattice->upperbounds = (mpz_t*)calloc(lgs_cols, sizeof(mpz_t));
        for (i = 0; i < lgs_cols; i++) {
            mpz_init_set_si(lattice->upperbounds[i], 1);
        }
        // Copy the upper bounds from the LGS and determine upperbounds_max,
        // which is the lcm of the non-zero upper bounds
        for (i = 0; i < lattice->num_boundedvars; i++) {
            mpz_set(lattice->upperbounds[i], LGS->upperbounds[i]);
            if (mpz_sgn(lattice->upperbounds[i]) != 0) {
                mpz_lcm(lattice->upperbounds_max, lattice->upperbounds_max, lattice->upperbounds[i]);
            }
        }
        if (mpz_cmp_si(lattice->upperbounds_max, 1) > 0) {
            lattice->is_zero_one = false;
        }
        fprintf(stderr, "upper bounds found. Max=");
        mpz_out_str(stderr, 10, lattice->upperbounds_max);
        fprintf(stderr, "\n");
    }
}

void handle_preselection(lgs_t *LGS, lattice_t *lattice) {
    int i;
    // Handle preselected columns
    if (LGS->original_cols != NULL) {
        lattice->num_cols_org = LGS->num_original_cols;
    } else {
        lattice->num_cols_org = LGS->num_cols;
    }

    lattice->original_cols = (int*)calloc(LGS->num_original_cols, sizeof(int));
    if (LGS->original_cols != NULL) {
        for (i = 0; i < lattice->num_cols_org; i++) {
            lattice->original_cols[i] = LGS->original_cols[i];
        }
    } else {
        for (i = 0; i < lattice->num_cols_org; i++) {
            lattice->original_cols[i] = 1;
        }
        fprintf(stderr, "No preselected columns \n");
    }
    fflush(stderr);

}

void init_diagonal_part(lgs_t *LGS, lattice_t *lattice) {
    int j;
    mpz_t upfac;

    mpz_init(upfac);

    // Append the other (diagonal) parts of lattice
    for (j = lattice->lgs_rows; j < lattice->num_rows; j++) {
        mpz_mul_si(
            lattice->basis[j - lattice->lgs_rows][j],
            lattice->max_norm,
            lattice->denom
        );
        mpz_mul_si(
            lattice->basis[lattice->num_cols - 1][j],
            lattice->max_norm,
            lattice->nom
        );
    }
    mpz_set(lattice->basis[lattice->lgs_cols + lattice->free_RHS][lattice->num_rows - 1], lattice->max_norm);

    if (lattice->free_RHS) {
        mpz_set_si(lattice->basis[lattice->lgs_cols][lattice->num_rows - 2], 1);
        mpz_set_si(lattice->basis[lattice->lgs_cols + 1][lattice->num_rows - 2], 0);
    }
    mpz_set(
        lattice->basis[lattice->lgs_cols + lattice->free_RHS][lattice->num_rows - 1], lattice->max_norm
    );

    // Multiply the diagonal entries and
    // the last columns to ensure the upper bounds on the variables
    mpz_init_set(lattice->max_norm_initial, lattice->max_norm);
    mpz_init_set_si(lattice->max_up, 1);

    if (!lattice->is_zero_one){
        for (j = 0; j < lattice->num_boundedvars; j++) {
            if (mpz_sgn(lattice->upperbounds[j]) != 0) {
                mpz_divexact(upfac, lattice->upperbounds_max, lattice->upperbounds[j]);
            } else {
                mpz_mul(upfac, lattice->upperbounds_max, lattice->upperbounds_max);
                mpz_mul_si(upfac, upfac, 10000);
            }
            mult_by(lattice->basis,
                    j, j + lattice->lgs_rows,
                    upfac);
            mult_by(lattice->basis,
                    lattice->lgs_cols + lattice->free_RHS, j + lattice->lgs_rows,
                    lattice->upperbounds_max);
        }
        mpz_set(lattice->max_up, lattice->upperbounds_max);
        mpz_mul(lattice->max_norm, lattice->max_norm, lattice->max_up);

        if (lattice->free_RHS) {
            mult_by(lattice->basis,
                    lattice->lgs_cols, lattice->num_rows - 2,
                    lattice->max_up);
        }

        mult_by(lattice->basis,
                lattice->lgs_cols + lattice->free_RHS, lattice->num_rows - 1,
                lattice->max_up);
    }
    mpz_clear(upfac);
}

void lgs_to_lattice(lgs_t *LGS, lattice_t *lattice) {
    int i, j;
    int lgs_rows = LGS->num_rows;
    int lgs_cols = LGS->num_cols;
    mpz_t factor;

    // Set the lattice dimensions
    lattice->num_rows = lgs_rows + lgs_cols + 1;
    lattice->num_cols = lgs_cols + 1;
    lattice->num_boundedvars = LGS->num_boundedvars;
    lattice->lgs_rows = lgs_rows;
    lattice->lgs_cols = lgs_cols;
    lattice->lgs_rank = LGS->rank;
    lattice->num_rows_org = lattice->num_rows;
    lattice->num_cols_org = lattice->num_cols;

    if (lattice->free_RHS) {
        lattice->num_rows++;
        lattice->num_cols++;
        fprintf(stderr,"The RHS is free !\n");
    } else {
        fprintf(stderr,"The RHS is fixed !\n");
    }

    alloc_basis(lattice);
    lattice->work_on_long = false;

    // Copy the linear system to the basis.
    for (j = 0; j < lgs_rows; j++) {
        for (i = 0; i < lgs_cols; i++) {
            mpz_set(lattice->basis[i][j], LGS->matrix[j][i]);
        }
        mpz_set(lattice->basis[lgs_cols][j], LGS->rhs[j]);
    }

    // Determine lcm of upper bounds on the variables
    handle_upperbounds(LGS, lattice);

    //
    // Multiply upper part of the lattice basis by
    // matrix_factor * upperbounds_max
    //
    mpz_init(factor);
    mpz_mul(factor, lattice->matrix_factor, lattice->upperbounds_max);
    for (j = 0; j < lgs_rows; j++) {
        for (i = 0; i < lgs_cols; i++) {
            mpz_mul(lattice->basis[i][j], lattice->basis[i][j], factor);
        }
        mpz_mul(lattice->basis[lgs_cols][j], lattice->basis[lgs_cols][j], factor);
    }
    mpz_clear(factor);

    // Store removed variables to include them in solution vector
    handle_preselection(LGS, lattice);

    lattice->nom = 1;
    lattice->denom = 2;
    init_diagonal_part(LGS, lattice);

    alloc_decomp(lattice);
}

DOUBLE** alloc_double_matrix(size_t cols, size_t rows) {
    int i;
    size_t rows_aligned, rows_size;

    // **m is a list of pointers m[0], ..., m[cols-1]
    // m[0] is a 1-dim array containing the whole matrix
    // m[1], ..., m[cols-1] point to the start
    // of each column.
    // This allows to use blas_ddot2 on non-contiguous arrays
    DOUBLE **m = (DOUBLE**)aligned_alloc(ALIGN_SIZE, DO_ALIGN(cols * sizeof(DOUBLE*)));

    // Changes here have to be done in cutlattice(), too
    rows_size = DO_ALIGN(rows * sizeof(DOUBLE));
    rows_aligned = rows_size / sizeof(DOUBLE);
    m[0] = (DOUBLE*)aligned_alloc(ALIGN_SIZE, cols * rows_size);

    for (i = 0; i < cols * rows_aligned; i++) {
        m[0][i] = 0.0;
    }
    for (i = 1; i < cols; i++) {
        m[i] = (DOUBLE*)(m[0] + i * rows_aligned);
    }
    return m;
}

int alloc_decomp(lattice_t *lattice) {
    int j, m;
    size_t len;
    int cols = lattice->num_cols;
    int rows = lattice->num_rows;
    decomp_t *d = &(lattice->decomp);

    if ((rows < 1) || (cols < 1)) return 0;

    len = DO_ALIGN(cols * sizeof(DOUBLE));
    d->c = (DOUBLE*)aligned_alloc(ALIGN_SIZE, len);
    d->N = (DOUBLE*)aligned_alloc(ALIGN_SIZE, len);
    for (j = 0; j < cols; j++) { d->c[j] = 0.0; }
    for (j = 0; j < cols; j++) { d->N[j] = 0.0; }

    // Use contiguous memory for BLAS
    // Attention: the columns of these matrices have to be adapted in cutlattice
    d->mu = alloc_double_matrix(cols, rows);
    m = (rows > cols) ? rows : cols;
    d->bd = alloc_double_matrix(m, rows);

    // R, H and h_beta are only pointers to already existing arrays
    d->R = d->mu;
    d->H = d->bd;
    d->h_beta = d->c;

    return 1;
}

int free_decomp(decomp_t d) {
    free(d.mu[0]);  // This is enough, see alloc_matrix()
    free(d.bd[0]);
    free(d.mu);
    free(d.bd);

    free(d.N);
    free(d.c);

    return 1;
}

double log_potential(DOUBLE **R, int s, int z) {
    double d = 0.0;
    int i;

    for (i = 0; i < s; i++) {
        d += log(fabs(R[i][i])) * (s - i);
    }
    d *= 0.5;
    return d;
}

double orthogonality_defect(lattice_t *lattice, DOUBLE **R, int s, int z) {
    double defect = 0.0;
    int i;

    for (i = 0; i < s; i++)
        defect += log(dot_mpz(lattice->basis[i], lattice->basis[i], lattice->num_rows)) - log(R[i][i]);

    defect *= 0.5;
    return defect;
}

void handle_signals(lattice_t *lattice, DOUBLE **R) {
    if (PRINT_REQUIRED && R != NULL) {
        //print_lattice(lattice, stderr);
        print_lattice_stat(lattice, R);
        PRINT_REQUIRED = 0;
    }
    if (DUMP_REQUIRED) {
        dump_lattice(lattice);
        DUMP_REQUIRED = 0;
    }
}

void print_lattice_stat(lattice_t *lattice, DOUBLE **R) {
    int i;
    fprintf(stderr, "--------------------------------\n");
    for (i = 0; i < lattice->num_cols; i++) {
        fprintf(stderr, "%0.3lf ", R[i][i] * R[i][i]);
    }
    fprintf(stderr,"\n logD=%0.3lf\n", log_potential(R, lattice->num_cols, lattice->num_rows));
    fprintf(stderr,"--------------------------------\n");
    fflush(stderr);
}

void shufflelattice(lattice_t *lattice) {
    mpz_t *swap;
    int i, j, r;
    unsigned int s;

    #if false
        s = (unsigned)(time(0))*getpid();
    #else
        s = 1300964772;
    #endif
    fprintf(stderr, "Seed=%u\n",s);
    srand(s);

    for (j = 0; j < 10000; j++) {
        for (i = lattice->num_cols - 1; i > 0; i--) {
            r = rand() % i;
            swap = lattice->basis[r];
            lattice->basis[r] = lattice->basis[i];
            lattice->basis[i] = swap;
        }
    }
    return;
}

void copy_lattice_to_long(lattice_t *lattice) {
    int i, j;

    for (i = 0; i < lattice->num_cols; ++i) {
        for (j = 0; j < lattice->num_rows; ++j) {
            lattice->basis_long[i][j] = mpz_get_si(lattice->basis[i][j]);
        }
    }
}

void copy_lattice_to_mpz(lattice_t *lattice) {
    int i, j;

    for (i = 0; i < lattice->num_cols; ++i) {
        for (j = 0; j < lattice->num_rows; ++j) {
            mpz_set_si(lattice->basis[i][j], lattice->basis_long[i][j]);
        }
    }
}

void print_gsa(DOUBLE **R, int n, int start) {
    int i, m;
    int res;
    // DOUBLE b1;
    FILE* f = fopen("gsa.tmp", "w");
    FILE* f2 = fopen("gsa1.tmp", "w");

    m = n; //4 * n / 5;
    DOUBLE sx = 0.0;
    DOUBLE sx2 = 0.0;
    DOUBLE sy = 0.0;
    DOUBLE sxy = 0.0;
    for (i = 0; i < m; i++) {
        sx += i;
        sy += log2(R[i][i] * R[i][i]);
        // if (isnan(sy)) {
        //     fprintf(stderr, "%d: %lf %lf\n", i, R[i][i], log(R[i][i] * R[i][i]));
        //     break;
        // }
        sxy += i * log2(R[i][i] * R[i][i]);
        sx2 += i * i;
    }
    double a = (m * sxy - sx * sy) / (m * sx2 - sx * sx);
    double b = (sy - a * sx) / m;
    fprintf(stderr, "%d: %lf %lf\n", start, a, b);

    // b1 = log2(R[0][0] * R[0][0]);
    for (i = 0; i < n; i++) {
        // fprintf(f, "%d %lf\n", i, b1 - log(R[i][i] * R[i][i]));
        // fprintf(f2, "%d %lf\n", i, b1 - log(R[n-1][n-1] * R[n-1][n-1]));
        fprintf(f, "%d %lf\n", i, log2(R[i][i] * R[i][i]));
        fprintf(f2, "%d %lf\n", i, a * i + b);
    }
    fclose(f);
    fclose(f2);
    res = system("mv gsa.tmp gsa.out");
    if (res != 0) {
        fprintf(stderr, "mv gsa.tmp gsa.out failed.\n");
    }
    res = system("mv gsa1.tmp gsa1.out");
    if (res != 0) {
        fprintf(stderr, "mv gsa1.tmp gsa1.out failed!\n");
    }

}
