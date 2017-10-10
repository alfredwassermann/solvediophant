#include <signal.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <gmp.h>

#include "const.h"
#include "lgs.h"
#include "datastruct.h"
#include "lattice.h"
#include "lll.h"
#include "bkz.h"
#include "dio2.h"

#if defined(USEBLAS)
    #define BLAS 1
#else
    #define BLAS 0
#endif

#if BLAS
    #include "common.h"
    #include "cblas.h"
#endif

/**
 * global variables
 */
//mpz_t max_norm_initial;
//mpz_t max_up;
mpz_t dummy;

int MAX_DUAL_BOUNDS = 512;

mpz_t lastlines_factor;

long num_solutions;
int level_max;

mpz_t soltest_u;
mpz_t soltest_s;
mpz_t soltest_upfac;

static FILE* fp;

long diophant(lgs_t *LGS, lattice_t *lattice, FILE* solfile, int restart, char *restart_filename) {

    int i, j;
    int block_size;
    DOUBLE lD, lDnew;
    coeff_t *swap_vec;

    /**
     * Initialize some globals
     */
    mpz_init(lastlines_factor);
    mpz_init(soltest_u);
    mpz_init(soltest_s);
    mpz_init_set_ui(soltest_upfac, 1);

    #if BLAS
        fprintf(stderr, "Use OpenBLAS\n");
        //openblas_set_num_threads(8);
    #endif

    if (!preprocess(LGS)) {
        fprintf(stderr, "Total number of solutions: 0\n\n");
        return 0;
    }
    lgs_to_lattice(LGS, lattice);

    /**
     * open solution file
     */
    fp = solfile;
    if (lattice->LLL_params.silent) fprintf(fp, "SILENT\n");
    fflush(fp);

    #if 0
    printf("After scaling\n");
    print_lattice(lattice);
    #endif

    #if 1 // Do reduction
    /**
     * permute lattice columns
     */
    swap_vec = lattice->basis[lattice->num_cols-1];
    for (i = lattice->num_cols - 1; i > 0; i--)
        lattice->basis[i] = lattice->basis[i - 1];
    lattice->basis[0] = swap_vec;

    #if 0
    printf("After permute\n");
    print_lattice(lattice);
    #endif
    //shufflelattice(lattice);
    /**
     * first reduction
     */
    mpz_set_ui(lastlines_factor, 1);
    fprintf(stderr, "\n"); fflush(stderr);
    if (!restart) {
        lll(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_LOW, CLASSIC_LLL);

        #if 0
            printf("After first reduction\n");
            print_lattice(lattice);
        #endif
        /**
         * cut the lattice
         */
        if (cutlattice(lattice)) {
            fprintf(stderr, "First reduction successful\n"); fflush(stderr);
        } else {
            fprintf(stderr, "First reduction not successful\n"); fflush(stderr);
            return 0;
        }
        #if 0
            printf("After cutting\n");
            print_lattice(lattice);
        #endif

        #if 1
            shufflelattice(lattice);
            /**
             * second reduction
             */
            mpz_set_ui(lastlines_factor, 1);
            lll(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_MED, POT_LLL);
            lll(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGH, POT_LLL);
            fprintf(stderr, "Second reduction successful\n"); fflush(stderr);
        #endif


    } else {
        load_lattice(lattice, restart_filename);
        //print_lattice(lattice);
    }

    #if 0
        printf("After second reduction\n");
        print_lattice(lattice);
    #endif

    #if 1  // Third reduction
    /**
     * scale last rows
     */
    mpz_set(lastlines_factor, lattice->LLL_params.scalelastlinefactor);
    for (i = 0; i < lattice->num_cols; i++) {
        mpz_mul(lattice->basis[i][lattice->num_rows].c,
                lattice->basis[i][lattice->num_rows].c, lastlines_factor);
    }
    if (lattice->free_RHS) {
        for (i = 0; i < lattice->num_cols; i++) {
            mpz_mul(lattice->basis[i][lattice->num_rows - 1].c,
                    lattice->basis[i][lattice->num_rows - 1].c, lastlines_factor);
        }
    }

    /**
     * third reduction
     */
    //fprintf(stderr, "\n"); fflush(stderr);

    if (lattice->LLL_params.iterate) {
        iteratedlll(lattice, lattice->num_cols, lattice->num_rows, lattice->LLL_params.iterate_no, LLLCONST_HIGH, POT_LLL);
    } else {
        //shufflelattice(lattice);

        for (block_size = 4; block_size <= lattice->LLL_params.bkz.beta; block_size += 4) {
            lD = lDnew;
            //shufflelattice(lattice);
            lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGHER,
                        block_size, lattice->LLL_params.bkz.p,
                        solutiontest, solutiontest_long);
            fprintf(stderr, "BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);
        }
    }
    fprintf(stderr, "Third reduction successful\n"); fflush(stderr);

    dump_lattice(lattice);

    /* undo scaling of last rows */
    for (i = 0; i < lattice->num_cols; i++) {
        mpz_divexact(lattice->basis[i][lattice->num_rows].c, lattice->basis[i][lattice->num_rows].c, lastlines_factor);
    }
    if (lattice->free_RHS) {
        for (i = 0; i < lattice->num_cols; i++) {
            mpz_divexact(lattice->basis[i][lattice->num_rows-1].c,
                lattice->basis[i][lattice->num_rows-1].c,
                lastlines_factor);
        }
    }
    #endif // Third reduction
    #else
    read_NTL_lattice();
    #endif // Do reduction

    if (lattice->LLL_params.print_ntl) {
        fprintf(stderr, "Print lattice for NTL and exit\n");
        print_NTL_lattice(lattice);   /* Version for the NTL output */
        return 0;
    }

    /**
     * explicit enumeration
     */
    fprintf(stderr, "\n"); fflush(stderr);
    num_solutions = explicit_enumeration(lattice, lattice->num_cols, lattice->num_rows);

    /**
     * close solution file;
     */
    if (lattice->LLL_params.silent)
        print_num_solutions(num_solutions);

    /**
     * free multiprecision memory
     */
    for (j = 0; j < lattice->num_cols/* + ADDITIONAL_COLS*/; j++) {
        for (i = 0; i <= lattice->num_rows; i++) {
            mpz_clear(lattice->basis[j][i].c);
        }
        free(lattice->basis[j]);
        free(lattice->basis_long[j]);
    }
    free(lattice->basis);
    free(lattice->basis_long);

    for (i = 0; i <= lattice->num_rows; i++) {
        mpz_clear(lattice->swap[i].c);
    }
    free(lattice->swap);
    mpz_clear(lattice->matrix_factor);
    mpz_clear(lattice->max_norm);

    mpz_clear(lastlines_factor);
    mpz_clear(soltest_u);
    mpz_clear(soltest_s);
    mpz_clear(soltest_upfac);

    //mpz_clear(max_norm_initial);
    //mpz_clear(max_up);
    //mpz_clear(upperbounds_max);

    // if (lattice->upperbounds != NULL) {
    //     for (i = 0; i < lattice->lgs_cols; i++) mpz_clear(lattice->upperbounds[i]);
    //     free(lattice->upperbounds);
    // }

    return num_solutions;
}

/**
 * Basic subroutines
 */
void print_num_solutions(long num_solutions) {
    fprintf(fp, "%ld solutions\n", num_solutions);
    fflush(fp);
}

int cutlattice(lattice_t *lattice) {
    int j, i, flag;
    int m;

    /**
     * delete unnecessary columns
     */
    j=0;
    do {
        if (lattice->basis[j][0].p > lattice->lgs_rows) {
            j++;
        } else {
            for (i = j + 1; i < lattice->num_cols; i++) {
                lattice->basis[i-1] = lattice->basis[i];
            }
            lattice->num_cols--;
        }
    } while (j < lattice->num_cols);


    /**
     * test for right hand side columns
     */
    flag = 0;
    for (i = 0; i < lattice->num_cols; i++)
        if (mpz_sgn(get_entry(lattice->basis, i, lattice->num_rows-1)) != 0) {
            flag = 1;
            break;
        }

    if (flag == 0) {
        fprintf(stderr, "Nonhomogenous solution not possible.\n"); fflush(stderr);
        exit(2);

        return 0;  /* Just for the compiler */
    }

    /* Now the rows are deleted. */
    for (j = 0; j < lattice->num_cols; j++)  {
       if (lattice->num_boundedvars == 0) {
            for (i = lattice->lgs_rows; i < lattice->num_rows; i++)
                put_to(lattice->basis, j, i - lattice->lgs_rows,
                    get_entry(lattice->basis, j, i));
        } else {
            for (i = lattice->lgs_rows; i < lattice->lgs_rows + lattice->num_boundedvars; i++)
                put_to(lattice->basis, j, i-lattice->lgs_rows,
                    get_entry(lattice->basis, j, i));
            for (i = lattice->lgs_rows + lattice->lgs_cols; i < lattice->num_rows; i++)
                put_to(lattice->basis, j,
                    i - lattice->lgs_rows - lattice->lgs_cols + lattice->num_boundedvars,
                    get_entry(lattice->basis, j, i));
        }
    }
    lattice->num_rows -= lattice->lgs_rows;
    lattice->num_rows -= (lattice->lgs_cols - lattice->num_boundedvars);

    // New start positions for the vectors in mu and bd.
    // This is necessary for BLAS access to it.
    for (i = 1; i < lattice->num_cols; i++) {
        lattice->decomp.mu[i] = (DOUBLE*)(lattice->decomp.mu[0] + i * lattice->num_rows);
    }
    m = (lattice->num_rows > lattice->num_cols) ? lattice->num_rows : lattice->num_cols;
    for (i = 1; i < m; i++) {
        lattice->decomp.bd[i] = (DOUBLE*)(lattice->decomp.bd[0] + i * lattice->num_rows);
    }

    for (j = 0; j < lattice->num_cols; j++) {
        coeffinit(lattice->basis[j],lattice->num_rows);
    }

    return 1;
}

int solutiontest(lattice_t *lattice, int position) {
    int i,j;
    int low, up;
    int end;

    /* test the last two rows */
    up = lattice->num_rows - 1 - lattice->free_RHS;

    if (mpz_cmpabs(get_entry(lattice->basis, position, lattice->num_rows - 1), lattice->max_norm) !=0 ) return 0;
    if (mpz_sgn(get_entry(lattice->basis, position, up)) ==0 ) return 0;

    /* test, if column is a solution */
    low = 0;
    if (lattice->num_cols == lattice->lgs_cols + 1 + lattice->free_RHS) {
        for (i = 0; i < lattice->lgs_rows; i++)
            if (mpz_sgn(get_entry(lattice->basis, position, i))!=0) return 0;
        low = lattice->lgs_rows;
    }

    if (lattice->is_zero_one) {
        for (i = low; i < up; i++) {
            if (mpz_cmpabs(get_entry(lattice->basis, position, i), lattice->max_norm)!=0) return 0;
        }
    } else {
        for (i=low;i<up;i++) {
            if (mpz_cmpabs(get_entry(lattice->basis, position,i), lattice->max_norm)>0) return 0;
        }
    }


    //mpz_set_si(upfac, 1);
    //mpz_divexact(soltest_s, get_entry(lattice->basis, position,
    //    lattice->num_rows - 1), lattice->LLL_params.scalelastlinefactor);
    mpz_set(soltest_s, get_entry(lattice->basis, position, lattice->num_rows - 1));

    /* write a solution with blanks */
    i = low;

    end = (lattice->cut_after == -1) ? lattice->no_original_cols : lattice->cut_after;

    for (j = 0; j < end; j++) {
        if (lattice->original_cols[j] == 0) {
            mpz_set_si(soltest_u, 0);
        } else {
            if (!lattice->is_zero_one) {
                if (mpz_cmp_si(lattice->upperbounds[i-low], 0) != 0) {
                    mpz_divexact(soltest_upfac, lattice->upperbounds_max, lattice->upperbounds[i - low]);
                } else {
                    mpz_set(soltest_upfac,lattice-> upperbounds_max);
                }
            }
            mpz_set(soltest_u, get_entry(lattice->basis, position, i));
            mpz_sub(soltest_u, soltest_u, soltest_s);
            mpz_divexact(soltest_u, soltest_u, lattice->max_norm_initial);
            mpz_divexact(soltest_u, soltest_u, soltest_upfac);
            mpz_divexact_ui(soltest_u, soltest_u, lattice->denom);
            mpz_abs(soltest_u, soltest_u);
            i++;
        }
        mpz_out_str(stderr, 10, soltest_u);
        fprintf(stderr, " ");
        if (lattice->LLL_params.stop_after_solutions == 1) {
            mpz_out_str(fp, 10, soltest_u);
            fprintf(fp," ");
        }
    }
    if (lattice->free_RHS) {
        mpz_divexact(soltest_u, get_entry(lattice->basis, position, up), lattice->max_up);
        // mpz_divexact(soltest_u, soltest_u, lattice->LLL_params.scalelastlinefactor);
        mpz_abs(soltest_u, soltest_u);
        fprintf(stderr, " L = ");
        mpz_out_str(stderr, 10, soltest_u);
    }
    fprintf(stderr, " !!\n");
    fflush(stderr);

    /* test if one solution is enough */
    if (lattice->LLL_params.stop_after_solutions == 1) {
        fprintf(fp,"\n");

        fprintf(stderr, "Stopped in phase 1 after finding a random solution\n");
        exit(8);
    }

    return 1;
}

int solutiontest_long(lattice_t *lattice, int position) {
    int i, j;
    int low, up;
    int end;

    #if 0
    int is_good = TRUE;
    for (j = 0; j < lattice->num_rows; ++j) {
        if (labs(lattice->basis_long[position][j]) != 1) {
            is_good = FALSE;
            break;
        }
    }
    if (is_good) {
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SOLUTION\n");
    }
    #endif

    for (j = 0; j < lattice->num_rows; ++j) {
        mpz_set_si(lattice->basis[position][j+1].c, lattice->basis_long[position][j]);
    }
    coeffinit(lattice->basis[position], lattice->num_rows);

    up = lattice->num_rows - 1 - lattice->free_RHS;

    /* test the last two rows */
    if (mpz_cmpabs(get_entry(lattice->basis, position, lattice->num_rows-1), lattice->max_norm) !=0 ) return 0;
    if (mpz_sgn(get_entry(lattice->basis, position, up)) ==0 ) return 0;

    /* test, if column is a solution */
    low = 0;
    if (lattice->num_cols == lattice->lgs_cols + 1 + lattice->free_RHS) {
        for (i = 0; i < lattice->lgs_rows; i++)
            if (mpz_sgn(get_entry(lattice->basis, position, i))!=0) return 0;
        low = lattice->lgs_rows;
    }

    if (lattice->is_zero_one) {
        for (i = low; i < up; i++) {
            if (mpz_cmpabs(get_entry(lattice->basis, position, i), lattice->max_norm) !=0) return 0;
        }
    } else {
        for (i = low; i < up; i++) {
            if (mpz_cmpabs(get_entry(lattice->basis, position, i), lattice->max_norm) > 0) return 0;
        }
    }

    // mpz_set_si(upfac, 1);
    // mpz_divexact(soltest_s, get_entry(lattice->basis, position, lattice->num_rows-1), lattice->LLL_params.scalelastlinefactor);
    mpz_set(soltest_s, get_entry(lattice->basis, position, lattice->num_rows - 1));

    /* write a solution with blanks */
    i = low;
    end = (lattice->cut_after == -1) ? lattice->no_original_cols : lattice->cut_after;
    for (j = 0; j < end; j++) {
        if (lattice->original_cols[j] == 0) {
            mpz_set_si(soltest_u,0);
        } else {
            if (lattice->is_zero_one) {
                if (mpz_cmp_si(lattice->upperbounds[i - low], 0) != 0) {
                    mpz_divexact(soltest_upfac, lattice->upperbounds_max, lattice->upperbounds[i-low]);
                } else {
                    mpz_set(soltest_upfac, lattice->upperbounds_max);
                }
            }
            mpz_set(soltest_u,get_entry(lattice->basis, position, i));
            mpz_sub(soltest_u, soltest_u, soltest_s);
            mpz_divexact(soltest_u, soltest_u, lattice->max_norm_initial);
            mpz_divexact(soltest_u, soltest_u, soltest_upfac);
            mpz_divexact_ui(soltest_u, soltest_u, lattice->denom);
            mpz_abs(soltest_u, soltest_u);
            i++;
        }
        mpz_out_str(stderr,10,soltest_u);
        fprintf(stderr, " ");
        if (lattice->LLL_params.stop_after_solutions == 1) {
            mpz_out_str(fp,10,soltest_u);
            fprintf(fp," ");
        }
    }
    if (lattice->free_RHS) {
        mpz_divexact(soltest_u, get_entry(lattice->basis, position, up), lattice->max_up);
        // mpz_divexact(soltest_u, soltest_u, lattice->LLL_params.scalelastlinefactor);
        mpz_abs(soltest_u, soltest_u);
        fprintf(stderr, " L = ");
        mpz_out_str(stderr, 10,soltest_u);
    }
    fprintf(stderr, " ||\n");
    fflush(stderr);

    /* test if one solution is enough */
    if (lattice->LLL_params.stop_after_solutions == 1) {
        fprintf(fp,"\n");

        fprintf(stderr, "Stopped in phase 1 after finding a random solution\n");
        exit(8);
    }

    return 1;
}


/**
 * LLL variants
 */
void lll(lattice_t *lattice, int s, int z, DOUBLE quality, int reduction_type) {
    DOUBLE **R = lattice->decomp.R;
    DOUBLE *beta = lattice->decomp.c;
    //DOUBLE *N = lattice->decomp.N;
    DOUBLE **H = lattice->decomp.H;
    int r, bit_size;

    //decomp_alloc(lattice); //&R, &beta, &N, &H, s, z);
    bit_size = get_bit_size(lattice);
    r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, solutiontest);
    //lllfree(R, beta, N, H, s);

    return;
}

DOUBLE iteratedlll(lattice_t *lattice, int s, int z, int no_iterates, DOUBLE quality, int reduction_type) {
    DOUBLE **R = lattice->decomp.R;
    DOUBLE *beta = lattice->decomp.c;
    //DOUBLE *N = lattice->decomp.N;
    DOUBLE **H = lattice->decomp.H;
    int r, i, j, runs;
    int bit_size;
    coeff_t *swapvl;
    long *swap;
    DOUBLE lD;

    //decomp_alloc(&R, &beta, &N, &H, s, z);

    bit_size = get_bit_size(lattice);

    if (bit_size < 32) {
        copy_lattice_to_long(lattice);
        r = lllH_long(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, solutiontest_long);
    } else {
        r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, solutiontest);
    }
    lD = log_potential(R, s, z);
    fprintf(stderr, "   log(D)= %f\n", lD);
    fflush(stderr);

    for (runs = 1; runs < no_iterates; runs++) {
        #if 0
        for (j = s - 1; j > 0; j--) {
            for (l = j - 1; l >= 0; l--) {
                /*|if (N[l] < N[j]) {|*/    /* $<$ sorts 'in descending order.' */
                if (N[l] > N[j]) {    /* $>$ sorts 'in ascending order.' */
                    swapvl = b[l];
                    for (i = l + 1; i <= j; i++) b[i-1] = b[i];
                    b[j] = swapvl;
                }
            }
        }
        #endif
        // Shuffle
        for (j = 0; j < 100; j++) {
            for (i = s - 1; i > 0; i--) {
                r = rand() % i;
                swapvl = lattice->basis[i];
                lattice->basis[i] = lattice->basis[r];
                lattice->basis[r] = swapvl;

                swap = lattice->basis_long[i];
                lattice->basis_long[i] = lattice->basis_long[r];
                lattice->basis_long[r] = swap;

            }
        }
        if (bit_size < 32) {
            r = lllH_long(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, solutiontest_long);
        } else {
            r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, solutiontest);
        }
        lD = log_potential(R, s, z);
        fprintf(stderr, "%d: log(D)= %f\n", runs, lD);
        fflush(stdout);
    }
    if (bit_size < 32) {
        copy_lattice_to_mpz(lattice);
    }
    //lllfree(R, beta, N, H, s);

    return lD;
}

DOUBLE block_reduce(lattice_t *lattice, int s, int z, int block_size, DOUBLE quality, int reduction_type) {
    DOUBLE **R = lattice->decomp.R;
    DOUBLE *beta = lattice->decomp.c;
    //DOUBLE *N = lattice->decomp.N;
    DOUBLE **H = lattice->decomp.H;

    DOUBLE lD;
    int start = 0, up, size, bit_size;
    coeff_t **basis_org;

    //decomp_alloc(&R, &beta, &N, &H, s, z);
    bit_size = get_bit_size(lattice);

    //r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size);
    #if 1
    while (start < s) {
        fprintf(stderr, "Block reduce %d\n", start);
        up = start + block_size;
        up = (up > s) ? s : up;

        basis_org = lattice->basis;
        lattice->basis = &(lattice->basis[start]);
        size = (start + block_size > up) ? up - start : block_size;
        lllH(lattice, R, beta, H, 0, 0, size, z, quality, reduction_type, bit_size, solutiontest);
        lattice->basis = basis_org;
        start += block_size;
    }
    #else
    while (start < s) {
        fprintf(stderr, "Block reduce %d\n", start);
        up = start + block_size;
        up = (up > s) ? s : up;

        lllH(lattice, R, beta, H, start, start, up, z, quality, reduction_type, bit_size);
        start += block_size;
    }
    #endif
    //print_lattice(lattice);

    lD = log_potential(R, s, z);
    fprintf(stderr, "   log(D)= %f\n", lD);
    fflush(stderr);
    //lllfree(R, beta, N, H, s);

    return lD;
}

/**
 * Exhaustive enumeration
*/

/**
 * Globals for enumeration
 */
/*|mpz_t *upb,*lowb;|*/
long dual_bound_success;
DOUBLE dum1, dum2;

long only_zeros_no, only_zeros_success, hoelder_no, hoelder_success;
long hoelder2_success;
long cs_success;

typedef struct {
    DOUBLE diff;
    DOUBLE cs;
    //DOUBLE l1;
    DOUBLE y;
    int num;
    int pos;
    DOUBLE us;
    DOUBLE* w;
} enum_node_t;

typedef struct {
    int pos;
    int num;
    int is_leave_count;
    enum_node_t* nodes;
} enum_level_t;

typedef struct {
    long loops;

    long *delta;
    long *d;
    long *eta;
    long *v;

    DOUBLE *y;
    DOUBLE *cs;
    DOUBLE *us;

    DOUBLE *coeff;
    DOUBLE **w;
} zigzag_t;

void allocateZigzag(zigzag_t *zigzag, int columns, int rows, int level) {
    int i, l,
        c = columns + 1;

    zigzag->us = (DOUBLE*)calloc(c, sizeof(DOUBLE));
    zigzag->cs = (DOUBLE*)calloc(c, sizeof(DOUBLE));
    zigzag->y = (DOUBLE*)calloc(c, sizeof(DOUBLE));
    zigzag->delta = (long*)calloc(c, sizeof(long));
    zigzag->d = (long*)calloc(c, sizeof(long));
    zigzag->eta = (long*)calloc(c, sizeof(long));
    zigzag->v = (long*)calloc(c, sizeof(long));
    zigzag->coeff = (DOUBLE*)calloc(c, sizeof(DOUBLE));

    zigzag->w = (DOUBLE**)calloc(c, sizeof(DOUBLE*));
    for (i = 0; i < c; i++) {
        zigzag->w[i] = (DOUBLE*)calloc(rows, sizeof(DOUBLE));
    }

    /* initialize arrays */
    for (i = 0; i < c; i++) {
        zigzag->cs[i] = zigzag->y[i] = zigzag->us[i] = 0.0;
        zigzag->delta[i] = 0;
        zigzag->v[i] = 0;
        zigzag->eta[i] = zigzag->d[i] = 1;
        for (l = 0; l < rows; l++) {
            zigzag->w[i][l] = 0.0;
        }
    }

    zigzag->us[level] = 1;
    zigzag->v[level] = 1;

    zigzag->loops = 0;
}

void allocateEnum_data(enum_level_t** enum_data, DOUBLE **fipo, int columns, int rows) {
    int i, j;
    long k;

    (*enum_data) = (enum_level_t*)calloc(columns + 1, sizeof(enum_level_t));
    for (i = 0; i <= columns; i++) {
        k = 2 * ((long)(fipo[0][i]) + 1);
        k = (k > MAX_DUAL_BOUNDS) ? MAX_DUAL_BOUNDS : k;
        (*enum_data)[i].nodes = (enum_node_t*)calloc(k, sizeof(enum_node_t));
        for (j = 0; j < k; j++) {
            (*enum_data)[i].nodes[j].w = (DOUBLE*)calloc(rows, sizeof(DOUBLE));
        }
        (*enum_data)[i].num = 0;
        (*enum_data)[i].pos = 0;
        (*enum_data)[i].is_leave_count = 0;
    }
}

int enumLevel(enum_level_t* enum_data, zigzag_t* z, lattice_t* lattice,
                DOUBLE* bd_1norm, DOUBLE** fipo,
                int* first_nonzero_in_column, int* firstp,
                int level, int rows, int columns, int bit_size, int max_steps) {

    DOUBLE old_coeff, stepWidth;
    int i, isSideStep;
    int goto_back;
    int is_good;

    DOUBLE** bd = lattice->decomp.bd;
    DOUBLE* c = lattice->decomp.c;
    DOUBLE Fd = lattice->decomp.Fd;
    DOUBLE Fqeps = lattice->decomp.Fqeps;
    DOUBLE Fq = lattice->decomp.Fq;

    enum_level_t* ed = &(enum_data[level]);
    ed->num = 0;
    isSideStep = FALSE;
    do {
        /* increase loop counter */
        z->loops++;
        if ((lattice->LLL_params.stop_after_loops > 0) &&
            (lattice->LLL_params.stop_after_loops <= z->loops)) {
            return -1;
        }

        #if VERBOSE > -1
        if (z->loops % 100000000 ==0) {                 /*10000000*/
            fprintf(stderr, "%ld loops, solutions: %ld",
                z->loops, num_solutions);
            fprintf(stderr, ", dual bounds: %ld ", dual_bound_success);
            fprintf(stderr, "\n");
            fflush(stderr);
        }
        #endif

        handle_signals(lattice, NULL);

        goto_back = FALSE;
        is_good = TRUE;

        /* compute new |cs| */
        old_coeff = z->coeff[level];
        z->coeff[level] = z->us[level] + z->y[level];
        z->cs[level] = z->cs[level+1] + z->coeff[level] * z->coeff[level] * c[level];

        if (z->cs[level] >= Fd)  {
            goto_back = TRUE;
        } else if (fabs(z->coeff[level]) > bd_1norm[level]) {
            /* Use (1, -1, 0, ...) as values in Hoelder pruning */
            goto_back = TRUE;
            ++hoelder2_success;
        } else if (fabs(z->us[level]) > fipo[0][level] ||
                   fabs(z->us[level] + 1) > fipo[1][level]) {
            dual_bound_success++;
            is_good = FALSE;
        } else {
            if (isSideStep) {
                stepWidth = z->coeff[level] - old_coeff;
                compute_w2(z->w, bd, stepWidth, level, rows);
            } else {
                compute_w(z->w, bd, z->coeff[level], level, rows);
            }
            if (level > 0) {
                i = prune_only_zeros(lattice, z->w, level, rows, Fq, first_nonzero_in_column, firstp,
                                      bd, z->y, z->us, columns);
                if (i < 0) {
                    goto_back = TRUE;
                } else if (i > 0) {
                    is_good = FALSE;
                } else if (prune(z->w[level], z->cs[level], rows, Fqeps)) {
                    ++hoelder_success;
                    is_good = FALSE;
                    if (z->eta[level] == 1) {
                        goto_back = TRUE;
                    } else {
                        z->eta[level] = 1;
                        z->delta[level] *= -1;
                        if (z->delta[level] * z->d[level] >= 0) z->delta[level] += z->d[level];
                        z->us[level] = z->v[level] + z->delta[level];
                        isSideStep = TRUE;
                        continue;
                    }
                }
            }
        }

        if (goto_back) {
            return 0;
        } else if (is_good) {
            i = ed->num;

            ed->nodes[i].us = z->us[level];
            ed->nodes[i].cs = z->cs[level];
            #if 1
                memcpy(ed->nodes[i].w, z->w[level], sizeof(DOUBLE)*rows);
            #else
                for (j = 0; j < rows; j++) {
                    ed->nodes[i].w[j] = z->w[level][j];
                }
            #endif
            ed->nodes[i].y = z->y[level];

            ed->num++;
            if (max_steps >= 0 && ed->num >= max_steps) {
                return 0;
            }

            if (ed->num >= MAX_DUAL_BOUNDS) {
                fprintf(stderr, "enum_data too small! Exit\n");
                fflush(stderr);
                exit(1);
            }
            isSideStep = TRUE;
        } else {
            isSideStep = FALSE;
        }

        /*
            Side step: the next value in the same level is
            chosen.
        */
        if (z->eta[level] == 0) {
            z->delta[level] *= -1;
            if (z->delta[level] * z->d[level] >= 0) {
                z->delta[level] += z->d[level];
            }
        } else {
            z->delta[level] += z->d[level] * ((z->delta[level] * z->d[level] >= 0) ? 1: -1);
        }
        z->us[level] = z->v[level] + z->delta[level];
    } while (TRUE);

    return 0;
}

int dfs(enum_level_t* enum_data, zigzag_t* z, lattice_t* lattice,
                // DOUBLE** bd, DOUBLE* c,
                //DOUBLE Fd, DOUBLE Fqeps, DOUBLE Fq,
                DOUBLE* bd_1norm, DOUBLE** fipo,
                int* first_nonzero_in_column, int* firstp,
                int level, int rows, int columns, int bit_size, DOUBLE** mu_trans) {

    int j;
    //DOUBLE s;
    enum_level_t* ed = &(enum_data[level]);

    if (-1 == enumLevel(enum_data, z, lattice,
            //bd, c, Fd, Fqeps, Fq,
            bd_1norm, fipo,
            first_nonzero_in_column, firstp,
            level, rows, columns, bit_size, -1)) {

        return -1;
    }

    for (ed->pos = 0; ed->pos < ed->num; ed->pos++) {
        z->us[level] = ed->nodes[ed->pos].us;
        z->cs[level] = ed->nodes[ed->pos].cs;

        // z->w[level] = ed->nodes[ed->pos].w;
        for (j = 0; j < rows; j++) {
            z->w[level][j] = ed->nodes[ed->pos].w[j];
        }

        if (level == 0) {
            // Solution found
            if (final_test(z->w[0], rows, lattice->decomp.Fq, z->us, lattice, bit_size) == 1) {
                print_solution(lattice, z->w[level], rows, lattice->decomp.Fq, z->us, columns);

                for (j = columns - 1 ; FALSE && j >= 0; j--) {
                    //if (1 || z->us[j] != ROUND(-z->y[j])) {
                    if (enum_data[j].pos > 0) {
                        fprintf(stderr, "====== ");
                    }
                    fprintf(stderr, "%d: %d of %d:\n",
                        j, //z->us[j], ROUND(-z->y[j]),
                        enum_data[j].pos, enum_data[j].num
                    );
                }

                if (lattice->LLL_params.stop_after_solutions > 0 &&
                    lattice->LLL_params.stop_after_solutions <= num_solutions) {
                    return -1;
                }
            }
        } else {

            level--;

            z->delta[level] = z->eta[level] = 0;
            z->y[level] = compute_y(mu_trans, z->us, level, level_max);
            z->us[level] = z->v[level] = ROUND(-z->y[level]);
            z->d[level] = (z->v[level] > -z->y[level]) ? -1 : 1;

            if (-1 == dfs(enum_data, z, lattice,
                //bd, c, Fd, Fqeps, Fq,
                bd_1norm, fipo,
                first_nonzero_in_column, firstp,
                level, rows, columns, bit_size, mu_trans)) {
                return -1;
            }

            level++;
        }
    }

    #if TRUE
        level++;
        if (level > level_max) {
            // If we reach a new level_max, a side step ahs to be done
            // in order to initialise z->us
            level_max = level;
            z->delta[level] += z->d[level] * ((z->delta[level] * z->d[level] >= 0) ? 1: -1);
            z->us[level] = z->v[level] + z->delta[level];
        }
    #endif
    return level;
}

int lds(enum_level_t* enum_data, zigzag_t* z, lattice_t* lattice,
                DOUBLE* bd_1norm, DOUBLE** fipo,
                int* first_nonzero_in_column, int* firstp,
                int level, int rows, int columns, int bit_size, DOUBLE** mu_trans,
                int lds_k, int lds_l, int lds_threshold) {

    int j;
    //int result;
    int start, end, pos, do_left_branch_last, p;
    int next_lds_k;
    int height, max_height, count;
    int max_steps;
    enum_level_t* ed = &(enum_data[level]);
    // DOUBLE** bd,
    // DOUBLE* c,
    // DOUBLE Fd,
    // DOUBLE Fqeps,
    // DOUBLE Fq,

    max_steps = -1;
    if (level >= lds_threshold && lds_k == 0) {
        max_steps = 1;
    }

    if (-1 == enumLevel(enum_data, z, lattice,
            //bd, c, Fd, Fqeps, Fq,
            bd_1norm, fipo,
            first_nonzero_in_column, firstp,
            level, rows, columns, bit_size, max_steps)) {

        return -1;
    }

    if (ed->num == 0) {
        ed->is_leave_count++;
    }

    start = 1;
    do_left_branch_last = 1;
    if (level < lds_threshold) {
        // dfs branching
        start = 0;
        end = ed->num;
        do_left_branch_last = 0;
    } else {
        // lds branching
        if (level - lds_threshold < lds_k) {
            // depth <= k -> no left branch
            start = 1;
            do_left_branch_last = 0;
        }
        if (lds_k > 0) {
            // Take all nodes per level:
            end = (lds_k < ed->num) ? lds_k + 1 : ed->num;
            // Take only two nodes per level:
            //end = (lds_k < 2) ? lds_k + 1 : 2;
        } else {
            // left-branches only
            #if 0
                // BBS
                // lds_k == 0: start conventional backtracking
                start = 0;
                end = ed->num;
                do_left_branch_last = 0;
            #else
                end = 1;
            #endif
        }
    }

    // BBS
    count = 0;
    max_height = 0;
    // for (pos = start; pos < end && pos < ed->num; pos++) {
    for (pos = start; pos <= ed->num; pos++) {
        // Right branches first
        if (pos >= end &&
            !(do_left_branch_last && pos == ed->num)) {
                continue;
            }
        p = pos % ed->num;
        ed->pos = p;
        //--------------
        // BBS
        // if (lds_k == 0 && count > 0) {
        //     break;
        // }
        //--------------
        z->us[level] = ed->nodes[p].us;
        z->cs[level] = ed->nodes[p].cs;
        // z->w[level] = ed->nodes[p].w;
        for (j = 0; j < rows; j++) {
            z->w[level][j] = ed->nodes[p].w[j];
        }

        if (level == 0) {
            // Solution found
            if (final_test(z->w[0], rows, lattice->decomp.Fq, z->us, lattice, bit_size) == 1) {
                print_solution(lattice, z->w[level], rows, lattice->decomp.Fq, z->us, columns);

                for (j = columns - 1 ; j >= 0; j--) {
                    fprintf(stderr, "%d: %d of %d\t%0.0lf\t%d",
                        j,
                        enum_data[j].pos, enum_data[j].num,
                        z->us[j],
                        enum_data[j].is_leave_count
                    );
                    if (enum_data[j].pos > 0) {
                        fprintf(stderr, "\t*");
                    }
                    fprintf(stderr, "\n");
                    // for (i = 0; i <= enum_data[j].num - 1; i++) {
                    //     s = enum_data[j].nodes[i].y + enum_data[j].nodes[i].us;
                    //     fprintf(stderr, "\t%.0lf\t%lf\t%lf\t%lf\t coeff=%lf\n",
                    //         enum_data[j].nodes[i].us,
                    //         enum_data[j].nodes[i].cs,
                    //         enum_data[j].nodes[i].l1,
                    //         s * s * c[j],
                    //         s
                    //     );
                    // }
                    if (j == lds_threshold) {
                        fprintf(stderr, "-------------------------------------\n");
                    }
                }

                if (lattice->LLL_params.stop_after_solutions > 0 &&
                    lattice->LLL_params.stop_after_solutions <= num_solutions)
                    return -1;
            }
        } else {
            level--;

            z->delta[level] = z->eta[level] = 0;
            //enum_data[level].pos = 0;
            z->y[level] = compute_y(mu_trans, z->us, level, level_max);
            z->us[level] = z->v[level] = ROUND(-z->y[level]);
            z->d[level] = (z->v[level] > -z->y[level]) ? -1 : 1;
            
            next_lds_k = lds_k;
            if (level >= lds_threshold) {
                // we are in ILDS mode
                if (p == 0) {
                    // depth > k, left branch
                    next_lds_k = lds_k;
                } else if (lds_k > 0) {
                    next_lds_k = (lds_k > p) ? lds_k - p : 0;
                }
            }
            
            height = lds(enum_data, z, lattice,
                    //bd, c, Fd, Fqeps, Fq,
                    bd_1norm, fipo,
                    first_nonzero_in_column, firstp,
                    level, rows, columns, bit_size, mu_trans, next_lds_k, lds_l, lds_threshold);

                    if (height == -1) {
                return -1;
            }
            if (height + 1 > max_height) max_height = height + 1;
            if (height >= lds_l) {
                count++;
            }
            
            level++;
        }

    }

    level++;
    if (level >= columns) {
        // We are done, let's leave the loop.
        //break;
        return 1;
    } else if (level > level_max) {
        level_max = level;
    }
    return max_height;
}

void init_dualbounds(lattice_t *lattice, DOUBLE ***fipo) {
    DOUBLE **muinv;
    DOUBLE entry;
    DOUBLE norm_1, norm_2;
    DOUBLE norm_1_1, norm_1_2;
    DOUBLE norm_2_1, norm_2_2;

    int i, j, l;
    int cols = lattice->num_cols;
    int rows = lattice->num_rows;

    (*fipo) = (DOUBLE**)calloc(cols + 1, sizeof(DOUBLE*));
    for (i = 0; i <= cols; i++) {
        (*fipo)[i] = (DOUBLE*)calloc(cols + 1, sizeof(DOUBLE));
    }

    muinv = (DOUBLE**)calloc(cols, sizeof(DOUBLE*));
    for(i = 0; i < cols; ++i) {
        muinv[i] = (DOUBLE*)calloc(rows, sizeof(DOUBLE));
    }

    /* determine inverse of mu */
    inverse(lattice->decomp.mu, muinv, cols);

    #if VERBOSE > -1
        fprintf(stderr, "Dual bounds:\n");
        fflush(stderr);
    #endif

    /* Symmetric Fincke-Pohst */
    for (i = 0; i < cols; i++) {
        norm_1 = norm_2 = 0.0;
        norm_1_1 = norm_1_2 = 0.0;
        norm_2_1 = norm_2_2 = 0.0;
        for (j = 0; j < rows; j++) {
            entry = 0.0;
            for (l = i; l < cols; l++) {
                entry += muinv[i][l] * lattice->decomp.bd[l][j] / lattice->decomp.c[l];
            }
            norm_2 += entry * entry;
            norm_1 += fabs(entry);
            #if TRUE
            for (l = cols - 1; l < cols; l++) {
                entry += muinv[cols - 1][l] * lattice->decomp.bd[l][j] / lattice->decomp.c[l];
            }
            norm_1_2 += entry * entry;
            norm_1_1 += fabs(entry);
            #endif
            #if TRUE
            for (l = cols - 2; l < cols; l++) {
                entry += muinv[cols - 2][l] * lattice->decomp.bd[l][j] / lattice->decomp.c[l];
            }
            norm_2_2 += entry * entry;
            norm_2_1 += fabs(entry);
            #endif

        }
        norm_2 = SQRT(norm_2 * lattice->decomp.Fd);
        norm_1 =  fabs(norm_1 * lattice->decomp.Fq) * (1.0 + EPSILON);
        (*fipo)[0][i] = (norm_1 < norm_2) ? norm_1 : norm_2;
        norm_1_2 = SQRT(norm_1_2 * lattice->decomp.Fd);
        norm_1_1 =  fabs(norm_1_1 * lattice->decomp.Fq) * (1.0 + EPSILON);
        (*fipo)[1][i] = (norm_1_1 < norm_1_2) ? norm_1_1 : norm_1_2;
        norm_2_2 = SQRT(norm_2_2 * lattice->decomp.Fd);
        norm_2_1 =  fabs(norm_2_1 * lattice->decomp.Fq) * (1.0 + EPSILON);
        (*fipo)[2][i] = (norm_2_1 < norm_2_2) ? norm_2_1 : norm_2_2;

        #if VERBOSE > -1
            fprintf(stderr, "%0.3lf ", (*fipo)[0][i]);
        #endif
    }

    // for (i = cols - 2; i >= 0; --i) {
    //     for (j = 0, tmp = 0.0; j < rows; j++) {
    //         dum1 = dual_basis[i][j] + dual_basis[i + 1][j];
    //         tmp += fabs(dum1);
    //     }
    //     dual_bound[i] = tmp * Fq * (1.0 + EPSILON);
    // }

    #if VERBOSE > -1
        fprintf(stderr, "\n\n");
        fflush(stderr);
    #endif
}

DOUBLE explicit_enumeration(lattice_t *lattice, int columns, int rows) {
    /* local variables for |explicit_enumeration() */
    /*|__attribute((aligned(16)))|*/

    int level;
    int i, j, l, k;
    int result;

    zigzag_t zigzag;
    enum_level_t* enum_data;

    DOUBLE *bd_1norm;
    int *first_nonzero, *first_nonzero_in_column, *firstp;
    int bit_size;

    DOUBLE **mu_trans;
    // DOUBLE *N, **mu, *c, **bd, **mu_trans;
    // DOUBLE **mu = lattice->decomp.R;
    // DOUBLE *c = lattice->decomp.c;
    // DOUBLE *N = lattice->decomp.N;
    // DOUBLE **bd = lattice->decomp.H;

    coeff_t *swap_vec;

    DOUBLE **fipo;
    /* test the size of the basis */
    fprintf(stderr, "Dimension of solution space (k): %d compared to columns-rank: %d\n",
                columns, lattice->lgs_cols - lattice->lgs_rank + 1 + lattice->free_RHS);
    fflush(stderr);

    if (columns < lattice->lgs_cols - lattice->lgs_rank + 1 + lattice->free_RHS) {
        fprintf(stderr,"LLL didn't succeed in computing a basis of the kernel.\n");
        fprintf(stderr,"Please increase c0 (the first parameter)!\n");
        return 0;
    }

    /* allocate the memory for enumeration */
    //decomp_alloc(lattice); //&mu, &c, &N, &bd, columns, rows);
    bd_1norm = (DOUBLE*)calloc(columns + 1, sizeof(DOUBLE));

    first_nonzero = (int*)calloc(rows, sizeof(int));
    first_nonzero_in_column = (int*)calloc(columns+rows+1, sizeof(int));
    if (first_nonzero_in_column == NULL) {
        return(0);
    }
    firstp = (int*)calloc(columns+1, sizeof(int));

    mu_trans = (DOUBLE**)calloc(columns+1, sizeof(DOUBLE*));
    for (i = 0; i <= columns; i++) {
        mu_trans[i]=(DOUBLE*)calloc(columns+1, sizeof(DOUBLE));
    }

    bit_size = get_bit_size(lattice);

    /* count nonzero entries in the last rows(s) */
    if (lattice->free_RHS) {
        i=0;
        for (j = columns - 1; j >= 0; j--) if (mpz_sgn(get_entry(lattice->basis, j, rows-2)) != 0)
            i++;
        fprintf(stderr, "Number of nonzero entries in the second last row: %d\n", i);
        fflush(stderr);
    }

    i = 0;
    for (j = columns - 1; j >= 0; j--) if (mpz_sgn(get_entry(lattice->basis, j, rows-1)) !=0 )
        i++;
    fprintf(stderr, "Number of nonzero entries in the last row: %d\n", i);
    fprintf(stderr, "Max bit size: %d\n", bit_size);
    fflush(stderr);

    // Move basis columns which have a nonzero entry in the last row to the end.
    // For lds this is mandatory.
    if (lattice->LLL_params.exhaustive_enum.lds == 1) {
        for (j = columns - 1; j > 0; j--) {
            for (l = j - 1; l >= 0;  l--) {
                if (mpz_cmpabs(get_entry(lattice->basis, l, rows - 1), get_entry(lattice->basis, j, rows - 1)) > 0) {
                    swap_vec = lattice->basis[l];
                    for (i = l + 1; i <= j; i++) lattice->basis[i - 1] = lattice->basis[i];
                    lattice->basis[j] = swap_vec;
                }
            }
        }
        //print_lattice(lattice);
    }

    /* set the simple pruning bounds */
    lattice->decomp.Fq = (DOUBLE)mpz_get_d(lattice->max_norm);
    lattice->decomp.Fd = (rows * lattice->decomp.Fq * lattice->decomp.Fq) * (1.0 + EPSILON);
    lattice->decomp.Fqeps = (1.0 + EPSILON) * lattice->decomp.Fq;        // Used in prune()
    #if VERBOSE > 0
        fprintf(stderr, "Fq: %f\n", (double)lattice->decomp.Fq);
        fprintf(stderr, "Fd: %f\n", (double)lattice->decomp.Fd);
        fflush(stderr);
    #endif

    /* orthogonalize the basis */
    #if GIVENS
        givens(lattice, columns, rows, lattice->decomp.mu, lattice->decomp.bd, lattice->decomp.c);
    #else
        gramschmidt(lattice, columns, rows, lattice->decomp.mu, lattice->decomp.bd, lattice->decomp.c);
    #endif

    /* compute $mu^\top$, the transpose of $mu$. */
    for (i = 0; i < columns; i++)
        for (j = 0; j < columns; j++)
            mu_trans[j][i] = lattice->decomp.mu[i][j];

    /* Compute 1-norm of orthogonal basis */
    for (i = 0; i <= columns; ++i) {
        bd_1norm[i] = 0.0;
        for (j = 0; j < rows; ++j) {
            bd_1norm[i] += fabs(lattice->decomp.bd[i][j]);
        }
        bd_1norm[i] *= lattice->decomp.Fqeps / lattice->decomp.c[i];
    }

    dual_bound_success = 0;
    init_dualbounds(lattice, &fipo);

    /* Remove trailing unnecessary columns.
     *
     * Contradiction to sorting columns, see above!
     * That means, columns whose corresponding Finke-Pohst bounds
     * are equal to 0 can be removed.
     * This is important for the Selfdual Bent Functions Problems
     */
    #if 1
    for (i = columns - 1; i >= 0; i--) {
        if (fipo[0][i] < 0.9) {
            fprintf(stderr, "DEL\n");
            columns--;
        } else {
            break;
        }
    }
    #endif

    /* New strategy */
    allocateEnum_data(&enum_data, fipo, columns, rows);

    /* initialize first-nonzero arrays */
    for (l = 0; l < rows; l++) {
        for (i = 0; i < columns; i++) if (mpz_sgn(get_entry(lattice->basis, i, l)) != 0) {
            first_nonzero[l] = i;
            break;
        }
    }

    fprintf(stderr, "First non-zero entries:\n");
    j = 0;
    for (l = 0; l < columns; l++) {
        firstp[l] = j;
        first_nonzero_in_column[j] = 0;
        j++;
        for (i = 0; i < rows; i++) {
            if (first_nonzero[i] == l) {
                first_nonzero_in_column[j] = i;
                first_nonzero_in_column[firstp[l]]++;
                j++;
            }
        }
        fprintf(stderr, "%d ", first_nonzero_in_column[firstp[l]]);
    }
    fprintf(stderr, ": %d\n", rows);
    firstp[columns] = j;
    first_nonzero_in_column[j] = 0;

    /* more initialization */
    level = first_nonzero[rows-1];
    if (level < 0) level = 0;

    level = first_nonzero[rows - 1];//columns - 1;
    level_max = level;

    allocateZigzag(&zigzag, columns, rows, level);

    only_zeros_no = only_zeros_success = 0;
    hoelder_no = hoelder_success = hoelder2_success = 0;
    cs_success = 0;

    fprintf(stderr, "Start enumeration at level=%d\n", level); fflush(stderr);
    /* the loop of the exhaustive enumeration */
    if (lattice->LLL_params.exhaustive_enum.lds == 1) {
        //for (i = 0; i <= columns / 2; i++) {
        //for (k = 0; k <= 8/*columns*/; k++) {
        for (k = 0; k < lattice->LLL_params.exhaustive_enum.lds_k_max; k++) {
            fprintf(stderr, "lds_k=%d\n", k); fflush(stderr);
            result = lds(enum_data, &zigzag, lattice,
                bd_1norm, fipo,
                first_nonzero_in_column, firstp,
                level, rows, columns, bit_size, mu_trans,
                k, 0, 0);

            if (result  == -1) {
                fprintf(stderr, "solution for lds_k=%d\n\n", k);
                break;
            }
        }
    } else {
        while (0 <= level && level < columns) {
            level = dfs(enum_data, &zigzag, lattice,
                //bd, c, Fd, Fqeps, Fq,
                bd_1norm, fipo,
                first_nonzero_in_column, firstp,
                level, rows, columns, bit_size, mu_trans);
        }
    }

    /* final output */
    fprintf(stderr, "Prune_cs: %ld\n", cs_success);
    fprintf(stderr, "Prune_only_zeros: %ld of %ld\n", only_zeros_success, only_zeros_no);
    fprintf(stderr, "Prune_hoelder: %ld of %ld\n", hoelder_success, hoelder_no);
    fprintf(stderr, "Prune_hoelder interval: %ld\n", hoelder2_success);
    fprintf(stderr, "Dual bounds: %ld\n", dual_bound_success);
    fprintf(stderr, "Loops: %ld\n", zigzag.loops);

    if ((lattice->LLL_params.stop_after_solutions <= num_solutions &&
         lattice->LLL_params.stop_after_solutions > 0) ||
        (lattice->LLL_params.stop_after_loops <= zigzag.loops &&
         lattice->LLL_params.stop_after_loops > 0 )) {
        fprintf(stderr, "Stopped after number of solutions: %ld\n", num_solutions);

        if (lattice->LLL_params.silent)
            print_num_solutions(num_solutions);
        if ((lattice->LLL_params.stop_after_loops <= zigzag.loops &&
            lattice->LLL_params.stop_after_loops > 0)) {
            exit(10);
        } else {
            exit(9);
        }
    } else {
        fprintf(stderr, "Total number of solutions: %ld\n", num_solutions);
    }
    fprintf(stderr, "\n");
    fflush(stdout);
    fflush(stderr);

    /* free allocated memory for enumeration */
    // free(us);
    // free(cs);
    // free(bd_1norm);
    // free(y);
    // free(delta);
    // free(d);
    // free(first_nonzero);
    // free(first_nonzero_in_column);
    // free(firstp);
    //
    // free(eta);
    // free(v);
    // for (l = 0; l <= columns; l++) free(w[l]);
    // free(w);
    // free(original_columns);

    // free(fipo);
    // for (l = 0; l < columns; l++) free(muinv[l]);
    // free(muinv);
    //
    // for (l = 0; l <= columns; l++) free(dual_basis[l]);
    // free(dual_basis);
    // free(dual_bound);

    // lllfree(mu, c, N, bd, columns);
    // for (l = 0; l < columns; l++) free(mu_trans[l]);
    // free(mu_trans);

    return 1;
}

DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max) {
    #if BLAS
        return cblas_ddot(level_max - level, &(us[level+1]), 1, &(mu_trans[level][level+1]), 1);
    #else
        int i;
        DOUBLE sum;
        i = level_max;
        sum = 0.0;
        while (i >= level + 1) {
            sum += mu_trans[level][i]*us[i];
            i--;
        }
        return sum;
    #endif
}

void compute_w2(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows) {
    #if BLAS
        cblas_daxpy(rows, alpha, bd[level], 1, w[level], 1);
    #else
        int i;
        for (i = rows - 1; i >= 0; --i) {
            w[level][i] += alpha * bd[level][i];
        }
    #endif

    return;
}

void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows) {
    #if BLAS
        cblas_dcopy(rows, w[level+1], 1, w[level], 1);
        if (fabs(alpha) > 1E-14) {
            cblas_daxpy(rows, alpha, bd[level], 1, w[level], 1);
        }
    #else
        int i;
        i = rows - 1;

        if (fabs(alpha) > 1E-14) {
            while (i >= 0) {
                w[level][i] = w[level+1][i] + alpha * bd[level][i];
                i--;
            }
        } else {
            while (i >= 0) {
                w[level][i] = w[level+1][i];
                i--;
            }
        }
    #endif

    return;
}

void gramschmidt(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c) {
    int i, l, j;
    DOUBLE sum;

    for (i = 0; i < columns; i++) {
        for (l = 0; l < rows; l++) bd[i][l] = (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l));
        for (j = 0; j < i; j++) {
            sum = 0.0;
            for (l = 0; l < rows; l++) sum += (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l)) * bd[j][l];
            mu[i][j] = sum / c[j];
            for (l = 0; l < rows; l++) bd[i][l] -= mu[i][j]*bd[j][l];
        }

        c[i] = scalarproductfp(bd[i], bd[i], rows);
        #if VERBOSE > 0
            fprintf(stderr, "%lf ",(double)c[i]);
        #endif
    }
    #if VERBOSE > 0
        fprintf(stderr, "\n\n");
        fflush(stderr);
    #endif
    return;
}

void givens(lattice_t *lattice, int columns, int rows, DOUBLE **mu,
            DOUBLE **bd, DOUBLE *c) {
    int i,l,j;
    int mm;
    DOUBLE d1,d2;
    DOUBLE gc,gs;
    DOUBLE t;


    /* The matrix |b| is copied to |mu|.
       |bd| is set to a $z\times z$ unity matrix.
    */
    for (i = 0; i < columns; i++) {
        for (l = 0; l < rows; l++) {
            mu[i][l] = (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l));
        }
    }

    for (i = 0; i < rows; i++) {
        for (l = 0; l < rows; l++) bd[i][l] = 0.0;
        bd[i][i] = 1.0;
    }

    for (j = 1; j < rows; j++) {    /* The Givens rotation */
        mm = (j < columns) ? j : columns;
        for (i = 0; i < mm; i++) {
            if (mu[i][j] == 0.0) {
                /* Nothing has to be done */
                gc = 1.0;
                gs = 0.0;
            } else {
                /* Stable computation of the
                   rotation coefficients.
                */
                if (fabs(mu[i][j]) >= fabs(mu[i][i])) {
                    t = mu[i][i] / mu[i][j];
                    gs = 1.0 / SQRT(1.0 + t*t);
                    gc = gs * t;
                } else {
                    t = mu[i][j] / mu[i][i];
                    gc = 1.0 / SQRT(1.0 + t*t);
                    gs = gc * t;
                }
                /* Rotation of |mu| */
                for (l = i; l < columns; l++) {
                    d1 = mu[l][i];
                    d2 = mu[l][j];
                    mu[l][i] =  gc*d1 + gs*d2;
                    mu[l][j] = -gs*d1 + gc*d2;
                }
                /* Rotation of the matrix $Q^t$ */
                for (l = 0; l < rows; l++) {
                    d1 = bd[i][l];
                    d2 = bd[j][l];
                    bd[i][l] =  gc*d1 + gs*d2;
                    bd[j][l] = -gs*d1 + gc*d2;
                }
            }
        }
    }

    /* Finally some scaling has to be done, since $Q$ is a orthonormal matrix */
    for (i = 0; i < columns; i++) {
        c[i] = mu[i][i] * mu[i][i];
        for (j = 0; j < rows; j++) {
            bd[i][j] *= mu[i][i];
        }
        for (j = columns - 1; j >= i; j--)
            mu[j][i] /= mu[i][i];

        #if VERBOSE > -1
            fprintf(stderr, "%6.3f ",(double)c[i]);
            if (i>0 && i%15==0) fprintf(stderr, "\n");
        #endif
    }
    #if VERBOSE > -1
        fprintf(stderr, "\n\n");
        fflush(stderr);
    #endif

    return;
}

void inverse(DOUBLE **mu, DOUBLE **muinv, int columns) {
    int i, j, k;
    DOUBLE sum;

    for (j = 0; j < columns; j++)
        for (i = j; i >= 0; i--) {
            sum = 0.0;
            for (k = i + 1; k < columns; k++)
                sum += mu[k][i]*muinv[k][j];
            if (i == j)
                muinv[i][j] = 1.0 - sum;
            else
                muinv[i][j] = -sum;
        }
    return;
}

/* There are several pruning methods.*/
int final_test(DOUBLE *v, int rows, DOUBLE Fq, DOUBLE *us, lattice_t *lattice, int bit_size) {
    register int i;
    register int k;

    i = rows - 1;
    do {
        if (fabs(v[i]) > Fq + 0.5 + EPSILON) {
            return 0;
        }
        i--;
    } while (i>=0);

    // If the involved numbers are too big,
    // an exact test is done.
    if (bit_size < 27) {
        return 1;
    }
    for (i = 0; i < rows; i++) {
        if (!lattice->is_zero_one) {
            if (mpz_cmp_si(lattice->upperbounds[i], 0) != 0) {
                mpz_divexact(soltest_upfac, lattice->upperbounds_max, lattice->upperbounds[i]);
            } else {
                mpz_set(soltest_upfac, lattice->upperbounds_max);
            }
        }

        mpz_set_si(soltest_u,0);
        for (k = 0; k < lattice->num_cols; k++) {
            if (ROUND(us[k]) > 0) {
                mpz_addmul_ui(soltest_u, get_entry(lattice->basis, k, i), ROUND(us[k]));
            } else {
                mpz_submul_ui(soltest_u, get_entry(lattice->basis, k, i), -ROUND(us[k]));
            }
        }

        mpz_sub(soltest_u, soltest_u, soltest_s);
        mpz_divexact(soltest_u, soltest_u, lattice->max_norm_initial);
        mpz_divexact(soltest_u, soltest_u, soltest_upfac);
        mpz_divexact_ui(soltest_u, soltest_u, lattice->denom);
        mpz_abs(soltest_u, soltest_u);
        if (!lattice->is_zero_one && (mpz_cmp_si(soltest_u, 0) < 0 ||
            mpz_cmp(soltest_u, lattice->upperbounds[i]) > 0) ) {
            //fprintf(stderr," rounding error -> this is not a solution!\n");
            return 0;
        }
    }

    return 1;
}

/* Pruning according to H\"olders inequality */
int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps) {
    #if BLAS
        if (cs < Fqeps * cblas_dasum(rows, w, 1)) {
            return 0;
        }
    #else
        register DOUBLE reseite;
        register int i;

        reseite = 0.0; /*|-cs/Fqeps;|*/ /* | * (1-eps) | */
        i = rows - 1;
        do {
            reseite += fabs(w[i]);
            i--;
        } while (i >= 0);
        if (cs < Fqeps * reseite) return 0;
    #endif

    return 1;
}

int prune_only_zeros(lattice_t *lattice, DOUBLE **w, int level, int rows, DOUBLE Fq,
                     int *first_nonzero_in_column, int *firstp,
                     DOUBLE **bd, DOUBLE *y, DOUBLE *us, int columns) {
    int i;
    int f;
    DOUBLE u1, u2;

    only_zeros_no++;
    for (i=0; i<first_nonzero_in_column[firstp[level]]; i++) {
        f = first_nonzero_in_column[firstp[level]+1+i];
        u1 = ( Fq-w[level+1][f])/bd[level][f] - y[level];
        u2 = (-Fq-w[level+1][f])/bd[level][f] - y[level];

        if (lattice->is_zero_one) {
            if (fabs(u1-round(u1))>EPSILON && fabs(u2-round(u2))>EPSILON) {
                only_zeros_success++;
                return -1;
            }

            if ( fabs(fabs(w[level][f])-Fq) > EPSILON ) {
                only_zeros_success++;
                return 1;
            }

        } else {  /* Not zero-one */

            /* Here we have to be very conservative */
            if (u2-u1 <= 1.0 + EPSILON &&
                    fabs(w[level][f]) < UINT32_MAX &&
                    fabs(w[level][f] - round(w[level][f])) > 0.001) {
                only_zeros_success++;
                return -1;
            }

            if (fabs(w[level][f]) > Fq * (1+EPSILON)) {
                return 1;
            }
        }
    }
    return 0;
}

int print_solution(lattice_t *lattice, DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns) {
    int i,j,k;
    int upper;
    int end;

    /* Test again, if the vector is really a solution */
    if (fabs(fabs(w[rows-1]) - Fq) > 0.5*Fq*EPSILON)  {
        return 0;
    }
    upper = rows - 1 - lattice->free_RHS;
    if (lattice->free_RHS && fabs(w[upper]) > Fq * (1 + EPSILON)) {
        return 0;
    }

    if (!lattice->LLL_params.silent) {
        mpz_set_si(soltest_upfac,1);
        mpz_set_si(soltest_s,0);
        for (k=0;k<columns;k++) {
            if (ROUND(us[k])>0) {
                mpz_addmul_ui(soltest_s,get_entry(lattice->basis, k, rows-1), ROUND(us[k]));
            } else {
            mpz_submul_ui(soltest_s,get_entry(lattice->basis, k,rows-1), -ROUND(us[k]));
            }
        }

        i = 0;
        end = (lattice->cut_after == -1) ? lattice->no_original_cols : lattice->cut_after;
        for (j = 0; j < end; j++) {
            if (lattice->original_cols[j] == 0) {
                mpz_set_si(soltest_u,0);
            } else {
                if (!lattice->is_zero_one) {
                    if (mpz_cmp_si(lattice->upperbounds[i],0)!=0) {
                        mpz_divexact(soltest_upfac, lattice->upperbounds_max, lattice->upperbounds[i]);
                    } else {
                        mpz_set(soltest_upfac, lattice->upperbounds_max);
                    }
                }
                mpz_set_si(soltest_u,0);
                for (k = 0; k < columns; k++) {
                    if (ROUND(us[k])>0) {
                        mpz_addmul_ui(soltest_u,get_entry(lattice->basis, k, i), ROUND(us[k]));
                    } else {
                        mpz_submul_ui(soltest_u,get_entry(lattice->basis, k, i), -ROUND(us[k]));
                    }
                }
                mpz_sub(soltest_u, soltest_u, soltest_s);
                mpz_divexact(soltest_u, soltest_u, lattice->max_norm_initial);
                mpz_divexact(soltest_u, soltest_u, soltest_upfac);
                mpz_divexact_ui(soltest_u, soltest_u, lattice->denom);
                mpz_abs(soltest_u, soltest_u);

                i++;
            }
            mpz_out_str(NULL, 10, soltest_u);
            fflush(stdout);
            mpz_out_str(fp, 10, soltest_u);

            /* Meanwhile, all solution vectors are written with separating blanks. */
            /*|if (!lattice->iszeroone) { }|*/
            printf(" ");
            fprintf(fp, " ");
        }
        if (lattice->free_RHS) {
            mpz_set_d(soltest_u, ROUND(w[i]));
            mpz_divexact(soltest_u, soltest_u, lattice->max_up);
            mpz_abs(soltest_u, soltest_u);
            printf(" L = ");
            mpz_out_str(NULL,10, soltest_u);
        }
        printf("\n");
        fflush(stdout);
        fprintf(fp, "\n");
        fflush(fp);
    }

    num_solutions++;
    if (num_solutions%10000==0) {
        printf("%ld\n", num_solutions);
        fflush(stdout);
    }

    return 1;
}

void stop_program_sig(int sig) {
    if (sig != SIGALRM)
       return;

    fprintf(stderr, "Stopped after SIGALRM, number of solutions: %ld\n", num_solutions);
    if (!SILENT)
        print_num_solutions(num_solutions);

    exit(11);
}

void print_lattice_sig(int sig) {
    if (sig != SIGUSR1)
       return;

    PRINT_REQUIRED = 1;
}

void dump_lattice_sig(int sig) {
    if (sig != SIGUSR2)
       return;

    DUMP_REQUIRED = 1;
}

void print_NTL_lattice(lattice_t *lattice) {
    int i,j;

    fprintf(stderr, "%d %d\n", lattice->num_cols, lattice->num_rows);
    //printf("%d\n",system_rows);
    printf("\n[");
    for (i = 0; i < lattice->num_cols; i++) {
        printf("[");
        for (j = 0; j < lattice->num_rows; j++) {
            mpz_out_str(NULL, 10, lattice->basis[i][j+1].c);
            printf(" ");
        }
        printf("]");
        printf("\n");
    }
    printf("]\n");
    fflush(stdout);

    printf("\n");
    //printf("%d ", lattice->num_cols - 1);
    mpz_out_str(NULL, 10, lattice->upperbounds_max);
    printf("\n\n[");
    for (i = 0; i < lattice->num_rows - 1; i++) {
        mpz_out_str(NULL, 10, lattice->upperbounds[i]);
        printf(" ");
    }
    printf("]\n");
    fflush(stdout);

    return;
}
