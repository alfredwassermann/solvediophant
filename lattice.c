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

#if defined(USEBLAS)
    #define BLAS 1
#else
    #define BLAS 0
#endif

#if BLAS
    #include "OpenBLASsub/common.h"
    #include "OpenBLASsub/cblas.h"
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
            mpz_out_str(NULL, 10, get_entry(lattice->basis, i, j));
            printf(" ");
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
            mpz_out_str(f, 10, get_entry(lattice->basis, i, j));
            fprintf(f, " ");
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
    fscanf(f, "%d%d\n", &(lattice->num_rows), &(lattice->num_cols));
    fscanf(f, "%d\n", &(lattice->cut_after));
    fscanf(f, "%d\n", &(lattice->free_RHS));

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

long gcd(long n1, long n2) {
    long a, b, c;

    if (n1 > n2) {
        a = n1;
        b = n2;
    } else {
        a = n2;
        b = n1;
    }

    while ((c = a % b) > 0) {
        a = b;
        b = c;
    }
    return b;
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

int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N,  DOUBLE ***bs, int s, int z) {
    int i, m;

    if ((z < 1) || (s < 1)) return 0;

    (*c) = (DOUBLE*)calloc(s, sizeof(DOUBLE));
    (*N) = (DOUBLE*)calloc(s, sizeof(DOUBLE));

    // Use contiguous memory for BLAS
    (*mu) = (DOUBLE**)calloc(s, sizeof(DOUBLE*));
    (*mu)[0] = (DOUBLE*)calloc(s * z, sizeof(DOUBLE));
    for (i = 1; i < s; i++) {
        (*mu)[i] = (DOUBLE*)((*mu)[0] + i * z); //
    }

    m = (z > s) ? z : s;
    (*bs) = (DOUBLE**)calloc(m,sizeof(DOUBLE*));
    (*bs)[0] = (DOUBLE*)calloc(z * m, sizeof(DOUBLE));
    for (i = 1; i < m; i++) {
        (*bs)[i] = (DOUBLE*)((*bs)[0] + i * z);
    }

    return 1;
}

int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s) {
    free(bs[0]);
    free(bs);

    free(mu[0]);
    free(mu);
    free(N);
    free(c);

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
    printf("--------------------------------\n");
    for (i = 0; i < lattice->num_cols; i++) {
        printf("%0.3lf ", R[i][i] * R[i][i]);
    }
    printf("\n logD=%0.3lf\n", log_potential(R, lattice->num_cols, lattice->num_rows));
    printf("--------------------------------\n");
    fflush(stdout);
}

void shufflelattice(lattice_t *lattice) {
    coeff_t *tmp;
    int i, j, r;
    unsigned int s;

    #if 0
        s = (unsigned)(time(0))*getpid();
    #else
        s = 1300964772;
    #endif
    fprintf(stderr, "Seed=%u\n",s);
    srand(s);

    for (j = 0; j < 100; j++) {
        for (i = lattice->num_cols - 1; i > 0; i--) {
            r = rand() % i;
            tmp = lattice->basis[r];
            lattice->basis[r] = lattice->basis[i];
            lattice->basis[i] = tmp;
        }
    }
    return;
}
