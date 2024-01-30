#include <signal.h>
#include <stdio.h>
#include <time.h>
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
void print_lattice(lattice_t *lattice) {
    int i, j;

    for (i = 0; i <= lattice->num_cols; i++) {
        for (j = 0; j < lattice->num_rows; j++) {
            if (lattice->work_on_long) {
                printf("%ld ", lattice->basis_long[i][j]);
            } else {
                mpz_out_str(NULL, 10, get_entry(lattice->basis, i, j));
                printf(" ");
            }
        }
        printf("\n");
    }
    printf("\n");
    fflush(stdout);
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
            res = mpz_inp_str(lattice->basis[i][j + 1].c, f, 10);
            if (res == 0) {
                incorrect_input_file();
            }
        }
        coeffinit(lattice->basis[i], lattice->num_rows);
    }
    fclose(f);
    return;
}

void coeffinit(coeff_t *v, int z) {
    short r = 0;
    short i;

    for (i = z; i >= 0; i--) {
        v[i].p = r;
        if (mpz_sgn(v[i].c) != 0) r = i;
    }

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
DOUBLE scalarproductlfp (coeff_t *v, coeff_t *w) {
    DOUBLE erg;
    long t1, t2;
    coeff_t *vv, *ww;

    erg = 0.0;
    t1 = v[0].p;
    t2 = w[0].p;
    if ((t1 == 0) || (t2 == 0))
        return 0;

    do {
        if (t2>t1) {
            t1 = v[t2-1].p;
            if (t2!=t1) {
                if (t1==0) break;
                t2 = w[t2].p;
                if (t2==0) break;
            }
            else goto gleich;
        } else if (t2<t1) {
            t2 = w[t1-1].p;
            if (t2!=t1) {
                if (t2==0) break;
                t1 = v[t1].p;
                if (t1==0) break;
            }
            else goto gleich;
        } else {
 gleich:    vv = &(v[t1]);
            ww = &(w[t2]);
            erg += (DOUBLE)mpz_get_d(vv->c) * (DOUBLE)mpz_get_d(ww->c);
            t1 = vv->p;
            if (t1==0) break;
            t2 = ww->p;
            if (t2==0) break;
        }
    } while (1);

    return (erg);
}

DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w , int n) {
    #if BLAS
        return cblas_ddot(n, v, 1, w, 1);
    #else
        DOUBLE r;
        int i;
        r = 0.0;
        for (i = n - 1; i >= 0; i--) r += v[i] * w[i];
        return r;
    #endif
}

void allocate_basis (lattice_t *lattice) {
    int i, j;

    // Allocate memory for the gmp basis and the long basis
    lattice->basis = (coeff_t**)calloc(lattice->num_cols + ADDITIONAL_COLS, sizeof(coeff_t*));
    lattice->basis_long = (long**)calloc(lattice->num_cols + ADDITIONAL_COLS, sizeof(long*));
    for (j = 0; j < lattice->num_cols + ADDITIONAL_COLS; j++) {
        lattice->basis_long[j] = (long*)calloc((unsigned int)lattice->num_rows, sizeof(long));
        lattice->basis[j] = (coeff_t*)calloc((unsigned int)lattice->num_rows + 1, sizeof(coeff_t));
        for (i = 0; i <= lattice->num_rows; i++) {
            mpz_init(lattice->basis[j][i].c);
        }
    }

    // Allocate memory for swap vector
    lattice->swap = (coeff_t*)calloc(lattice->num_rows + 1, sizeof(coeff_t));
    lattice->swap_long = (long*)calloc(lattice->num_rows + 1, sizeof(long));
    for (i = 0; i <= lattice->num_rows; i++) {
        mpz_init(lattice->swap[i].c);
    }

}

void handle_upperbounds(lgs_t *LGS, lattice_t *lattice) {
    int i;
    int lgs_cols = LGS->num_cols;

    // Handle upper bounds
    mpz_init_set_si(lattice->upperbounds_max, 1);

    lattice->is_zero_one = TRUE;
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
            lattice->is_zero_one = FALSE;
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
        lattice->no_original_cols = LGS->num_original_cols;
    } else {
        lattice->no_original_cols = LGS->num_cols;
    }

    lattice->original_cols = (int*)calloc(LGS->num_original_cols, sizeof(int));
    if (LGS->original_cols != NULL) {
        for (i = 0; i < lattice->no_original_cols; i++) {
            lattice->original_cols[i] = LGS->original_cols[i];
        }
    } else {
        for (i = 0; i < lattice->no_original_cols; i++) {
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
        mpz_mul_si(lattice->basis[j - lattice->lgs_rows][j + 1].c, lattice->max_norm, lattice->denom);
        mpz_mul_si(lattice->basis[lattice->num_cols - 1][j + 1].c, lattice->max_norm, lattice->nom);
    }
    mpz_set(lattice->basis[lattice->lgs_cols + lattice->free_RHS][lattice->num_rows].c, lattice->max_norm);

    if (lattice->free_RHS) {
        mpz_set_si(lattice->basis[lattice->lgs_cols][lattice->num_rows - 1].c, 1);
        mpz_set_si(lattice->basis[lattice->lgs_cols + 1][lattice->num_rows - 1].c, 0);
    }
    mpz_set(lattice->basis[lattice->lgs_cols + lattice->free_RHS][lattice->num_rows].c, lattice->max_norm);

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
            mult_by(lattice->basis, j, j + lattice->lgs_rows, upfac);
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
}

void lgs_to_lattice(lgs_t *LGS, lattice_t *lattice) {
    int i, j;
    int lgs_rows = LGS->num_rows;
    int lgs_cols = LGS->num_cols;

    // Set the lattice dimensions
    lattice->num_rows = lgs_rows + lgs_cols + 1;
    lattice->num_cols = lgs_cols + 1;
    lattice->num_boundedvars = LGS->num_boundedvars;
    lattice->lgs_rows = lgs_rows;
    lattice->lgs_cols = lgs_cols;
    lattice->lgs_rank = LGS->rank;

    if (lattice->free_RHS) {
        lattice->num_rows++;
        lattice->num_cols++;
        fprintf(stderr,"The RHS is free !\n");
    } else {
        fprintf(stderr,"The RHS is fixed !\n");
    }

    allocate_basis(lattice);
    lattice->work_on_long = FALSE;

    // Copy the linear system to the basis.
    // Thereby, multiply the entries by a large (enough) factor.
    for (j = 0; j < lgs_rows; j++) {
        for (i = 0; i < lgs_cols; i++) {
            mpz_mul(lattice->basis[i][j+1].c, LGS->matrix[j][i], lattice->matrix_factor);
        }
        mpz_mul(lattice->basis[lgs_cols][j+1].c, LGS->rhs[j], lattice->matrix_factor);
    }

    handle_upperbounds(LGS, lattice);
    handle_preselection(LGS, lattice);

    lattice->nom = 1;
    lattice->denom = 2;

    init_diagonal_part(LGS, lattice);

    // Init sparse structure
    for (i = 0; i < lattice->num_cols; i++) {
        coeffinit(lattice->basis[i], lattice->num_rows);
    }
    coeffinit(lattice->swap, lattice->num_rows);

    decomp_alloc(lattice);
}

int decomp_alloc(lattice_t *lattice) {
    int i, m;
    int cols = lattice->num_cols;
    int rows = lattice->num_rows;
    decomp_t *d = &(lattice->decomp);

    if ((rows < 1) || (cols < 1)) return 0;

    d->c = (DOUBLE*)calloc(cols, sizeof(DOUBLE));
    d->N = (DOUBLE*)calloc(cols, sizeof(DOUBLE));

    // Use contiguous memory for BLAS
    d->mu = (DOUBLE**)calloc(cols, sizeof(DOUBLE*));
    d->mu[0] = (DOUBLE*)calloc(cols * rows, sizeof(DOUBLE));
    for (i = 1; i < cols; i++) {
        d->mu[i] = (DOUBLE*)(d->mu[0] + i * rows);
    }

    m = (rows > cols) ? rows : cols;
    d->bd = (DOUBLE**)calloc(m, sizeof(DOUBLE*));
    d->bd[0] = (DOUBLE*)calloc(rows * m, sizeof(DOUBLE));
    for (i = 1; i < m; i++) {
        d->bd[i] = (DOUBLE*)(d->bd[0] + i * rows);
    }

    // R, H and h_beta are only pointers to already existing arrays
    d->R = d->mu;
    d->H = d->bd;
    d->h_beta = d->c;

    return 1;
}

int decomp_free(lattice_t *lattice) {
    decomp_t d = lattice->decomp;

    free(d.bd[0]);
    free(d.bd);

    free(d.mu[0]);
    free(d.mu);
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
        defect += log(scalarproductlfp(lattice->basis[i], lattice->basis[i])) - log(R[i][i]);

    defect *= 0.5;
    return defect;
}

void handle_signals(lattice_t *lattice, DOUBLE **R) {
    if (PRINT_REQUIRED && R != NULL) {
        //print_lattice(lattice);
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
    coeff_t *tmp;
    int i, j, r;
    unsigned int s;

    #if TRUE
        s = (unsigned)(time(0))*getpid();
    #else
        s = 1300964772;
    #endif
    fprintf(stderr, "Seed=%u\n",s);
    srand(s);

    for (j = 0; j < 10000; j++) {
        for (i = lattice->num_cols - 1; i > 0; i--) {
            r = rand() % i;
            tmp = lattice->basis[r];
            lattice->basis[r] = lattice->basis[i];
            lattice->basis[i] = tmp;
        }
    }
    return;
}

void copy_lattice_to_long(lattice_t *lattice) {
    int i, j;

    for (i = 0; i < lattice->num_cols; ++i) {
        for (j = 0; j < lattice->num_rows; ++j) {
            lattice->basis_long[i][j] = mpz_get_si(lattice->basis[i][j+1].c);
        }
    }
}

void copy_lattice_to_mpz(lattice_t *lattice) {
    int i, j;

    for (i = 0; i < lattice->num_cols; ++i) {
        for (j = 0; j < lattice->num_rows; ++j) {
            mpz_set_si(lattice->basis[i][j+1].c, lattice->basis_long[i][j]);
        }
        coeffinit(lattice->basis[i], lattice->num_rows);
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
