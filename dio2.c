#include <signal.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <gmp.h>
#include "dio2.h"
#include "const.h"

#if defined(USEBLAS)
    #define BLAS 1
#else
    #define BLAS 0
#endif

#if BLAS
    #include "OpenBLASsub/common.h"
    #include "OpenBLASsub/cblas.h"
#endif

/**
 * global variables
 */
mpz_t max_norm_initial;
mpz_t max_up;
mpz_t dummy;

long nom, denom;
mpz_t lastlines_factor;

int system_rows, system_columns;
//int lattice->num_rows, lattice->num_cols;

int free_RHS;
int iszeroone;
mpz_t *upperbounds;
mpz_t upperbounds_max;
mpz_t upfac;

int *original_columns;
int no_original_columns;
long nosolutions;
int nboundvars;

mpz_t soltest_u;
mpz_t soltest_s;
mpz_t soltest_upfac;

static FILE* fp;

/**
 * Inline functions
 */
#define ROUND(r) ceil(r-0.5)
#define put_to(mat, i, j, val) mpz_set(mat[i][j+1].c, val)
#define smult_lattice(mat, i, j, factor) mpz_mul(mat[i][j+1].c, mat[i][j+1].c, factor)
#define get_entry(mat, i, j) mat[i][j+1].c


long diophant(lgs_t *LGS, lattice_t *lattice, FILE* solfile, int restart, char *restart_filename) {

    int i,j;
    DOUBLE lD, lDnew;
    coeff_t *swap_vec;


    /**
     * Initialize some globals
     */
    mpz_init(lastlines_factor);
    mpz_init(upfac);

    mpz_init(soltest_u);
    mpz_init(soltest_s);
    mpz_init_set_ui(soltest_upfac, 1);

    nom = 1;
    denom = 2;

    system_rows = LGS->num_rows;
    system_columns = LGS->num_cols;
    nboundvars = LGS->num_boundedvars;

#if BLAS
    fprintf(stderr, "Use OpenBLAS\n");
    //openblas_set_num_threads(8);
#endif

    /**
     * set the lattice dimension;
     */
    lattice->num_rows = system_rows + system_columns + 1;
    lattice->num_cols = system_columns + 1;

    if (free_RHS) {
        lattice->num_rows++;
        lattice->num_cols++;
    } else {
        fprintf(stderr,"The RHS is fixed !\n");
        fflush(stderr);
    }

    /**
     * allocate memory
     */
    lattice->basis = (coeff_t**)calloc(lattice->num_cols + ADDITIONAL_COLS, sizeof(coeff_t*));
    for (j = 0; j < lattice->num_cols/* + ADDITIONAL_COLS*/; j++) {
        lattice->basis[j] = (coeff_t*)calloc(lattice->num_rows + 1, sizeof(coeff_t));
        for (i = 0; i <= lattice->num_rows; i++)
            mpz_init(lattice->basis[j][i].c);
    }
    // lattice->basis_s = (coeff_t**)calloc(lattice->num_cols + ADDITIONAL_COLS, sizeof(coeff_t*));
    // for (j = 0; j < lattice->num_cols + ADDITIONAL_COLS; j++) {
    //     lattice->basis_s[j] = (coeff_t*)calloc(lattice->num_rows + 1, sizeof(coeff_t));
    //     for (i = 0; i <= lattice->num_rows; i++)
    //         mpz_init(lattice->basis_s[j][i].c);
    // }

    lattice->swap = (coeff_t*)calloc(lattice->num_rows + 1, sizeof(coeff_t));
    for (i = 0; i <= lattice->num_rows; i++)
        mpz_init(lattice->swap[i].c);

    /**
     * read the system
     */
    for (j = 0; j < system_rows; j++) {
        for (i = 0; i < system_columns; i++) {
            mpz_mul(lattice->basis[i][j+1].c, LGS->matrix[j][i], lattice->matrix_factor);
        }
        mpz_mul(lattice->basis[system_columns][j+1].c, LGS->rhs[j], lattice->matrix_factor);
    }

    /**
     * handle upper bounds
     */
    mpz_init_set_si(upperbounds_max,1);
    iszeroone = 1;
    if (LGS->upperbounds == NULL) {
        fprintf(stderr, "No upper bounds: 0/1 variables are assumed \n"); fflush(stderr);
    } else {
        upperbounds = (mpz_t*)calloc(system_columns, sizeof(mpz_t));
        for (i = 0; i < system_columns; i++)
            mpz_init_set_si(upperbounds[i], 1);
        for (i = 0; i < nboundvars/*|system_columns|*/; i++) {
            mpz_set(upperbounds[i], LGS->upperbounds[i]);
            if (mpz_sgn(upperbounds[i]) != 0) {
                mpz_lcm(upperbounds_max, upperbounds_max, upperbounds[i]);
            }
        }
        if (mpz_cmp_si(upperbounds_max, 1) > 0)
            iszeroone = 0;

        fprintf(stderr,"upper bounds found. Max=");
        fflush(stderr);

        mpz_out_str(stderr, 10, upperbounds_max);
        fprintf(stderr, "\n");
        fflush(stderr);
    }

    /**
     * handle preselected columns
     */
    if (LGS->original_cols != NULL)
        no_original_columns = LGS->num_original_cols;
    else
        no_original_columns = LGS->num_cols;

    original_columns = (int*)calloc(LGS->num_original_cols, sizeof(int));

    if (LGS->original_cols != NULL)
        for (i = 0; i < no_original_columns; i++)
            original_columns[i] = LGS->original_cols[i];
    else {
        for (i = 0; i < no_original_columns; i++)
            original_columns[i] = 1;
        fprintf(stderr, "No preselected columns \n");
        fflush(stderr);
    }

    /**
     * append the other parts of lattice
     */
    for (j = system_rows; j < lattice->num_rows; j++) {
        mpz_mul_si(lattice->basis[j-system_rows][j+1].c, lattice->max_norm, denom);
        mpz_mul_si(lattice->basis[lattice->num_cols-1][j+1].c, lattice->max_norm, nom);
    }
    mpz_set(lattice->basis[system_columns+free_RHS][lattice->num_rows].c, lattice->max_norm);

    if (free_RHS) {
        mpz_set_si(lattice->basis[system_columns][lattice->num_rows-1].c, 1);
        mpz_set_si(lattice->basis[system_columns+1][lattice->num_rows-1].c, 0);
    }
    mpz_set(lattice->basis[system_columns+free_RHS][lattice->num_rows].c, lattice->max_norm);
    for (i = 0; i < lattice->num_cols; i++)
        coeffinit(lattice->basis[i], lattice->num_rows);
    coeffinit(lattice->swap, lattice->num_rows);

    /**
     * open solution file
     */
    fp = solfile;
    if (lattice->LLL_params.silent) fprintf(fp,"SILENT\n");
    fflush(fp);

#if 0
    printf("Before scaling\n");
    print_lattice(lattice);
#endif
    /**
     * scale lattice
     */
    mpz_init_set(max_norm_initial, lattice->max_norm);
    mpz_init_set_si(max_up, 1);
    if (!iszeroone){
        for (j=0;j<nboundvars/*|system_columns|*/;j++) {
            if (mpz_sgn(upperbounds[j])!=0) {
                mpz_divexact(upfac,upperbounds_max,upperbounds[j]);
            } else {
                mpz_mul(upfac,upperbounds_max,upperbounds_max);
                mpz_mul_si(upfac, upfac, 10000);
            }
            smult_lattice(lattice->basis, j, j + system_rows, upfac);
            smult_lattice(lattice->basis, system_columns+free_RHS, j + system_rows, upperbounds_max);
        }
        mpz_set(max_up, upperbounds_max);
        mpz_mul(lattice->max_norm, lattice->max_norm, max_up);
        if (free_RHS)
            smult_lattice(lattice->basis, system_columns, lattice->num_rows-2, max_up);

        smult_lattice(lattice->basis, system_columns+free_RHS, lattice->num_rows-1, max_up);
    }

#if 0
    printf("After scaling\n");
    print_lattice();
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
        #if 0
            block_reduce(lattice, lattice->num_cols, lattice->num_rows, 50, LLLCONST_HIGHER, DEEPINSERT_CONST);
            block_reduce(lattice, lattice->num_cols, lattice->num_rows, 100, LLLCONST_HIGHER, DEEPINSERT_CONST);
            block_reduce(lattice, lattice->num_cols, lattice->num_rows, 200, LLLCONST_HIGHER, DEEPINSERT_CONST);
            block_reduce(lattice, lattice->num_cols, lattice->num_rows, 400, LLLCONST_HIGHER, DEEPINSERT_CONST);
        #endif

        lll(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_LOW, -1);

        #if 0
            printf("After first reduction\n");
            print_lattice();
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
            print_lattice();
        #endif

        #if 1
            shufflelattice(lattice);
            /**
             * second reduction
             */
            mpz_set_ui(lastlines_factor, 1);
            lll(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_MED, 10);
            lll(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGH, 10);
            fprintf(stderr, "Second reduction successful\n"); fflush(stderr);
        #endif


    } else {
        load_lattice(lattice, restart_filename);
        //print_lattice(lattice);
    }


#if 0
    printf("After second reduction\n");
    print_lattice();
#endif

#if 1  // Third reduction
    /**
     * scale last rows
     */
    mpz_set(lastlines_factor, lattice->LLL_params.scalelastlinefactor);
    for (i = 0; i < lattice->num_cols; i++)
        mpz_mul(lattice->basis[i][lattice->num_rows].c, lattice->basis[i][lattice->num_rows].c, lastlines_factor);
    if (free_RHS)
        for (i = 0; i < lattice->num_cols; i++)
            mpz_mul(lattice->basis[i][lattice->num_rows-1].c, lattice->basis[i][lattice->num_rows-1].c, lastlines_factor);

#if 0
    for (i=0;i<lattice->num_cols;i++) {
        for (j=0;j<40;j++)
            mpz_mul_ui(lattice->basis[i][j+1].c,lattice->basis[i][j+1].c, 9);
    }
#endif

    /**
     * third reduction
     */
    fprintf(stderr, "\n"); fflush(stderr);

    if (lattice->LLL_params.iterate) {
        iteratedlll(lattice, lattice->num_cols, lattice->num_rows, lattice->LLL_params.iterate_no, LLLCONST_HIGH, DEEPINSERT_CONST);
    } else {
        //shufflelattice(lattice);

        i = 0;
        do {
            lD = lDnew;
            //shufflelattice(lattice);
            lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGHER,
                            lattice->LLL_params.bkz.beta, lattice->LLL_params.bkz.p);
            fprintf(stderr, "BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);
            fflush(stderr);
            i++;
        }
        while (i < 1 && fabs(lDnew - lD) > 0.01);
    }
    fprintf(stderr, "Third reduction successful\n"); fflush(stderr);

    dump_lattice(lattice);

    /* undo scaling of last rows */
    for (i = 0; i < lattice->num_cols; i++)
        mpz_divexact(lattice->basis[i][lattice->num_rows].c, lattice->basis[i][lattice->num_rows].c, lastlines_factor);
    if (free_RHS)
        for (i = 0; i < lattice->num_cols; i++)
            mpz_divexact(lattice->basis[i][lattice->num_rows-1].c,
                lattice->basis[i][lattice->num_rows-1].c,
                lastlines_factor);
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
    nosolutions = explicit_enumeration(lattice,lattice->num_cols,lattice->num_rows);

    /**
     * close solution file;
     */
    if (lattice->LLL_params.silent)
        print_num_solutions(nosolutions);

    /**
     * free multiprecision memory
     */
    for (j = 0; j < lattice->num_cols; j++) {
        for (i = 0; i <= lattice->num_rows; i++) {
            mpz_clear(lattice->basis[j][i].c);
        }
        free(lattice->basis[j]);
    }
    free(lattice->basis);
    for (i = 0; i <= lattice->num_rows; i++) {
        mpz_clear(lattice->swap[i].c);
    }
    free(lattice->swap);
    mpz_clear(lattice->matrix_factor);
    mpz_clear(lattice->max_norm);

    mpz_clear(lastlines_factor);
    mpz_clear(upfac);
    mpz_clear(max_norm_initial);
    mpz_clear(max_up);
    mpz_clear(soltest_u);
    mpz_clear(soltest_s);
    mpz_clear(soltest_upfac);
    mpz_clear(upperbounds_max);

    if (upperbounds != NULL) {
        for (i = 0; i < system_columns; i++) mpz_clear(upperbounds[i]);
        free(upperbounds);
    }

    return nosolutions;
}

/**
 * Basic subroutines
 */
void print_num_solutions(long num_solutions) {
    fprintf(fp, "%ld solutions\n", nosolutions);
    fflush(fp);
}

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
    for (i = 0; i < lattice->num_cols; i++) {
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

int cutlattice(lattice_t *lattice) {
    int j, i, flag;

    /**
     * delete unnecessary columns
     */
    j=0;
    do {
        if (lattice->basis[j][0].p > system_rows)
            j++;
        else {
            for (i = j + 1; i < lattice->num_cols; i++)
                lattice->basis[i-1] = lattice->basis[i];
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
       if (nboundvars == 0) {
            for (i = system_rows; i < lattice->num_rows; i++)
                put_to(lattice->basis, j, i-system_rows,
                    get_entry(lattice->basis, j, i));
        } else {
            for (i = system_rows; i < system_rows + nboundvars; i++)
                put_to(lattice->basis, j, i-system_rows,
                    get_entry(lattice->basis, j, i));
            for (i = system_rows+system_columns; i < lattice->num_rows; i++)
                put_to(lattice->basis, j, i-system_rows-system_columns+nboundvars,
                    get_entry(lattice->basis, j, i));
        }
    }
    lattice->num_rows -= system_rows;
    lattice->num_rows -= (system_columns-nboundvars);

    for (j = 0; j < lattice->num_cols; j++) coeffinit(lattice->basis[j],lattice->num_rows);

    return 1;
}

int solutiontest(lattice_t *lattice, int position) {
    int i,j;
    int low, up;
    int end;

    /* test the last two rows */
    if (mpz_cmpabs(get_entry(lattice->basis, position, lattice->num_rows-1), lattice->max_norm) !=0 ) return 0;
    if (mpz_sgn(get_entry(lattice->basis, position, lattice->num_rows-1-free_RHS)) ==0 ) return 0;

    /* test, if column is a solution */
    low = 0;
    up = lattice->num_rows-1-free_RHS;
    if (lattice->num_cols == system_columns + 1 + free_RHS) {
        for (i = 0; i < system_rows; i++)
            if (mpz_sgn(get_entry(lattice->basis, position, i))!=0) return 0;
        low = system_rows;
    }

    if (iszeroone) {
        for (i=low;i<up;i++) {
            if (mpz_cmpabs(get_entry(lattice->basis, position,i), lattice->max_norm)!=0) return 0;
        }
    } else {
        for (i=low;i<up;i++) {
            if (mpz_cmpabs(get_entry(lattice->basis, position,i), lattice->max_norm)>0) return 0;
        }
    }


    mpz_set_si(upfac,1);
    mpz_divexact(soltest_s, get_entry(lattice->basis, position, lattice->num_rows-1), lattice->LLL_params.scalelastlinefactor);

    /* write a solution with blanks */
    i = low;

    end = (lattice->cut_after == -1) ? no_original_columns : lattice->cut_after;

    for (j = 0; j < end; j++) {
        if (original_columns[j] == 0) {
            mpz_set_si(soltest_u,0);
        } else {
            if (!iszeroone) {
                if (mpz_cmp_si(upperbounds[i-low],0) != 0) {
                    mpz_divexact(soltest_upfac,upperbounds_max,upperbounds[i-low]);
                } else {
                    mpz_set(soltest_upfac,upperbounds_max);
                }
            }
            mpz_set(soltest_u,get_entry(lattice->basis, position, i));
            mpz_sub(soltest_u,soltest_u,soltest_s);
            mpz_divexact(soltest_u,soltest_u,max_norm_initial);
            mpz_divexact(soltest_u,soltest_u,soltest_upfac);
            mpz_divexact_ui(soltest_u,soltest_u,denom);
            mpz_abs(soltest_u,soltest_u);
            i++;
        }
        mpz_out_str(NULL,10,soltest_u);
        printf(" ");
        if (lattice->LLL_params.stop_after_solutions == 1) {
            mpz_out_str(fp,10,soltest_u);
            fprintf(fp," ");
        }
    }
    if (free_RHS) {
        mpz_divexact(soltest_u, get_entry(lattice->basis, position, up), max_up);
        mpz_divexact(soltest_u, soltest_u, lattice->LLL_params.scalelastlinefactor);
        mpz_abs(soltest_u,soltest_u);
        printf(" L = ");
        mpz_out_str(NULL,10,soltest_u);
    }
    printf("\n");
    fflush(stdout);

    /* test if one solution is enough */
    if (lattice->LLL_params.stop_after_solutions == 1) {
        fprintf(fp,"\n");

        fprintf(stderr, "Stopped in phase 1 after finding a random solution\n");
        exit(8);
    }

    return 1;
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
 *  Lattice basis reduction algorithms
 */
int lllHfp(lattice_t *lattice, DOUBLE **R, DOUBLE *beta, DOUBLE **H,
            int start, int low, int up, int z,
            DOUBLE delta, int deepinsert_blocksize,
            int bit_size) {

    coeff_t **b = lattice->basis;
    int i, j, k;
    DOUBLE norm;
    int mu_all_zero;

    DOUBLE theta, eta;

    DOUBLE mus;
    mpz_t musvl;
    mpz_t hv;

    DOUBLE r_new, r_act;
    DOUBLE pot, pot_max;
    int pot_idx;
    int insert_pos, lowest_pos;
    coeff_t *swapvl;

    #if VERBOSE > 0
        int counter = 0;
    #endif

    eta = ETACONST;
    if (bit_size > 100) {
        theta = 0.50;
        eta = 0.52;
    } else if (bit_size > 55) {
        theta = 0.05;
    } else if (bit_size > 30) {
        theta = 0.01;
    } else {
        theta = 0.0;
    }
    fprintf(stderr, "Start LLLHfp with deepinsert %d; max bits: %d\n",  deepinsert_blocksize, bit_size);

    mpz_init(musvl);
    mpz_init(hv);
    //mpz_init(sum_mu);

    /* Test for trivial cases. */
    if ((z <= 1) || (up <= 1)) {
        fprintf(stderr, "Wrong dimensions in LLLHfp\n");
        fflush(stderr);
        return(0);
    }

    k = (start >= low) ? start : low;
    lowest_pos = k;

    /* The main loop */
    while (k < up) {
        if (k < lowest_pos) {
            lowest_pos = k;
        }
        #if VERBOSE > 0
            if ((counter % 10000) == 0) {
                fprintf(stderr, "LLL_H: %d k:%d\n", counter, k);
                fflush(stderr);
            }
            counter++;
        #endif
        handle_signals(lattice);

        #if 0
        // Look ahead
        i = householder_column(b, R, H, beta, k, s, z, bit_size);
        if (i > k) {
            swapvl = b[i];
            for (j = i; j > k; --j) {
                b[j] = b[j - 1];
            }
            b[k] = swapvl;

            //fprintf(stderr, "GET %d from %d\n", k, i);
        }
        #endif


    start_tricol:
        /* Recompute column k of R */
        i = householder_column(b, R, H, beta, k, k + 1, z, bit_size);

        /* size reduction of $b_k$ */
        mu_all_zero = 1;
        for (j = k - 1; j >= low; j--) {
            if (fabs(R[k][j]) > eta * fabs(R[j][j]) + theta * fabs(R[k][k])) {
                mus = ROUND(R[k][j] / R[j][j]);
                mpz_set_d(musvl, mus);
                mu_all_zero = 0;

                /* set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$ */
                size_reduction(b, R, musvl, mus, k, j);
                solutiontest(lattice, k);
            } else {
            }
        }
        if (!mu_all_zero) {
            goto start_tricol;
        }

        if (0) {
            check_precision(b[k], R[k], z, k);
        }

        /*
            Before going to step 4 we test if $b_k$ is linear dependent.
            If we find a linear dependent vector $b_k$,
            we shift b_k to the last column of the
            matrix and restart lllHfp with s = s-1.
        */
        #if BLAS
            norm = cblas_dnrm2(k + 1, R[k], 1);
        #else
            for (j = 0, norm = 0.0; j <= k; ++j) {
                norm += R[k][j] * R[k][j];
            }
            norm = SQRT(norm);
        #endif
        if (norm < 0.5) {
            swapvl = b[k];
            for (i = k + 1; i < up; i++) {
                b[i-1] = b[i];
            }
            b[up - 1] = swapvl;
            up--;
            k = 0;
            fprintf(stderr, "Zero vector at %d\n", k);
            continue;
        }

        // If delta == 0, only size reduction is done
        if (delta == 0.0) {
            k++;
            continue;
        }

        /* fourth step: swap columns */
        #if 0
        // Standard LLL or deepinsert
        if (deepinsert_blocksize > 0) {
            i = low;
            #if BLAS
                r_new = cblas_ddot(k + 1, R[k], 1, R[k], 1);
            #else
                for (j = 0, r_new = 0.0; j <= k; ++j) {
                    r_new += R[k][j] * R[k][j];
                }
            #endif
        } else {
            i = (k > low) ? k - 1 : low;
            r_new = R[k][k] * R[k][k] + R[k][i] * R[k][i];
        }

        insert_pos = k;
        while (i < k) {
            r_act = delta * R[i][i] * R[i][i];
            //if (delta * R[i][i]*R[i][i] > rhs) {
            if (0.8 * r_act > r_new ||
                ((i < deepinsert_blocksize || k - i < deepinsert_blocksize) &&
                  r_act > r_new)) {
                insert_pos = i;
                break;
            }
            r_new -= R[k][i]*R[k][i];
            i++;
        }
        /*
        if (insert_pos < k) {
            fprintf(stderr, "swap %d to %d\n", k, insert_pos);
        }
        */
        #else
        // Pot-LLL
        if (deepinsert_blocksize > 0) {
            pot = pot_max = 0.0;
            pot_idx = k;
            for (i = k - 1; i >= low; --i) {
                for (j = k, r_new = 0.0; j >= i; --j) {
                    r_new += R[k][j] * R[k][j];
                }
                pot += log(r_new) - log(R[i][i] * R[i][i]);

                if (pot < log(delta)/* && pot < pot_max*/) {
                    pot_max = pot;
                    pot_idx = i;
                }
            }

            // if (pot_idx < k) {
            //     fprintf(stderr, "swap %d to %d: gain=%lf\n", k, pot_idx, pot_max);
            // }
            insert_pos = pot_idx;
        } else {
            insert_pos = k;
            i = (k > low) ? k - 1 : low;
            r_new = R[k][k] * R[k][k] + R[k][i] * R[k][i];
            r_act = delta * R[i][i] * R[i][i];
            if (r_act > r_new) {
                insert_pos = i;
            }
        }
        #endif

        if (insert_pos < k) {
            swapvl = b[k];
            for (j = k; j > insert_pos; --j) {
                b[j] = b[j - 1];
            }
            b[insert_pos] = swapvl;

            //fprintf(stderr, "INSERT %d at %d\n", k, insert_pos);
            k = insert_pos;
        } else {
            k++;
        }
    }
    mpz_clear(hv);
    mpz_clear(musvl);

    return lowest_pos;

}

int householder_column(coeff_t **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z, int bit_size) {
    int i, j;
    int l;
    DOUBLE zeta;
    DOUBLE w, w_beta;
    DOUBLE norm;
    DOUBLE eps = 0.0000000001;
    #if !BLAS
        DOUBLE x;
    #endif

    DOUBLE min_val = 0.0;
    DOUBLE min_idx = -1;

    for (l = k; l < s; l++) {
        for (j = 0; j < z; ++j) {
            R[k][j] = (DOUBLE)mpz_get_d(b[l][j+1].c);
        }

    #if BLAS
        // Compute R[k]
        for (i = 0; i < k; ++i) {
            w = cblas_ddot(z - i, &(R[k][i]), 1, &(H[i][i]), 1);
            w_beta = w * beta[i];
            cblas_daxpy(z - i, -w_beta, &(H[i][i]), 1, &(R[k][i]), 1);
        }
        // |R[k]|
        if (bit_size < 27) {
            norm = cblas_dnrm2(z - k, &(R[k][k]), 1);
        } else {
            j = cblas_idamax(z - k, &(R[k][k]), 1);
            zeta = fabs(R[k][k + j]);
            cblas_dscal(z - k, 1 / zeta, &(R[k][k]), 1);
            norm = zeta * cblas_dnrm2(z - k, &(R[k][k]), 1);
            cblas_dscal(z - k, zeta, &(R[k][k]), 1);
        }

        // H[k] = R[k] / |R[k]|
        cblas_dcopy(z - k, &(R[k][k]), 1, &(H[k][k]), 1);
        cblas_dscal(z - k, 1 / norm, &(H[k][k]), 1);

        H[k][k] += (R[k][k] >= -eps) ? 1.0 : -1.0;
        beta[k] = 1.0 / (1.0 + fabs(R[k][k]) / norm);

        // w = <R[k], H[k]>
        w = cblas_ddot(z - k, &(R[k][k]), 1, &(H[k][k]), 1);
        w_beta = w * beta[k];

        // R[k] -= w * beta * H[k]
        cblas_daxpy(z - k, -w_beta, &(H[k][k]), 1, &(R[k][k]), 1);
    #else
        // Compute R[k]
        for (i = 0; i < k; ++i) {
            for (j = i, w = 0.0; j < z; ++j) {
                w += R[k][j] * H[i][j];
            }
            w_beta = w * beta[i];
            for (j = i; j < z; ++j) {
                R[k][j] -= w_beta * H[i][j];
            }
        }

        // |R[k]|
        if (bit_size < 27) {
            for (j = k, norm = 0.0; j < z; ++j) {
                norm += x * x;
            }
        } else {
            // Use zeta for stability
            for (j = k, zeta = 0.0; j < z; ++j) {
                if (fabs(R[k][j]) > zeta) {
                    zeta = fabs(R[k][j]);
                }
            }
            for (j = k, norm = 0.0; j < z; ++j) {
                x = R[k][j] / zeta;
                norm += x * x;
            }
            norm = zeta * SQRT(norm);
        }

        // H[k] = R[k] / |R[k]|
        for (j = k; j < z; ++j) {
            H[k][j] = R[k][j] / norm;
        }

        H[k][k] += (R[k][k] >= -eps) ? 1.0 : -1.0;
        beta[k] = 1.0 / (1.0 + fabs(R[k][k]) / norm);

        // w = <R[k], H[k]>
        for (j = k, w = 0.0; j < z; ++j) {
            w += R[k][j] * H[k][j];
        }

        // R[k] -= w * beta * H[k]
        w_beta = w * beta[k];
        for (j = k; j < z; ++j) {
            R[k][j] -= w_beta * H[k][j];
        }
    #endif
        if (l == k || R[k][k] * R[k][k] < min_val) {
            min_val = R[k][k] * R[k][k];
            min_idx = l;
        }
    }
    return min_idx;
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

void check_precision(coeff_t *b, DOUBLE *R, int z, int k) {
    int j;
    mpz_t b_norm;
    DOUBLE r_norm;

    mpz_init(b_norm);
    for (j = 0, r_norm = 0.0; j < z; ++j) {
        mpz_addmul(b_norm, b[j+1].c, b[j+1].c);
    }
    for (j = 0, r_norm = 0.0; j <= k; ++j) {
        r_norm += R[j] * R[j];
    }
    if (fabs(mpz_get_d(b_norm) - r_norm) > 0.1) {
        fprintf(stderr, "precision check fails at %d: ", k);
        mpz_out_str(stderr, 10, b_norm);
        fprintf(stderr, " %lf\n", r_norm);
        fflush(stderr);
    }
    mpz_clear(b_norm);
}

void size_reduction(coeff_t **b, DOUBLE  **mu, mpz_t musvl, double mus, int k, int j) {
    int i, ii, iii;
    coeff_t *bb;

    switch (mpz_get_si(musvl)) {
    case 1:
        /* $\lceil\mu_{k,j}\rfloor = 1$ */
        i = b[j][0].p;
        while (i != 0) {
                bb = &(b[k][i]);
                mpz_sub(bb->c, bb->c, b[j][i].c);
                iii = bb->p;
                if ((b[k][i-1].p != i) && (mpz_sgn(bb->c) != 0))
                    for (ii = i - 1; (ii >= 0) && (b[k][ii].p == iii); ii--)
                        b[k][ii].p = i;
                else if (mpz_sgn(bb->c) == 0) {
                    for (ii = i - 1;  (ii >= 0) && (b[k][ii].p == i); ii--)
                        b[k][ii].p = iii;
                }
                i = b[j][i].p;
        }
    #if BLAS
        cblas_daxpy(j, -1.0, mu[j], 1, mu[k], 1);
    #else
        for (i = 0; i < j; i++) mu[k][i] -= mu[j][i];
    #endif

        break;

    case -1:
        /* $\lceil\mu_{k,j}\rfloor = -1$ */
        i = b[j][0].p;
        while (i != 0) {
                bb = &(b[k][i]);
                mpz_add(bb->c, bb->c, b[j][i].c);
                iii = bb->p;
                if ((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0))
                    for (ii = i - 1; (ii >= 0) && (b[k][ii].p == iii); ii--) b[k][ii].p = i;
                else if (mpz_sgn(bb->c)==0) {
                    for (ii = i - 1; (ii >= 0) && (b[k][ii].p == i); ii--) b[k][ii].p = iii;
                }
                i = b[j][i].p;
        }

    #if BLAS
        cblas_daxpy(j, 1.0, mu[j], 1, mu[k], 1);
    #else
        for(i = 0; i < j; i++) mu[k][i] += mu[j][i];
    #endif
        break;

    default:
        /* $\lceil\mu_{k,j}\rfloor \neq \pm 1$ */
        i = b[j][0].p;
        while (i != 0) {
                bb = &(b[k][i]);
                mpz_submul(bb->c, b[j][i].c, musvl);
                iii = bb->p;
                if ((b[k][i-1].p != i) && (mpz_sgn(bb->c) != 0))
                    for (ii = i - 1; (ii >= 0) && (b[k][ii].p ==  iii); ii--) b[k][ii].p = i;
                else if (mpz_sgn(bb->c) == 0) {
                    for (ii = i - 1; (ii >= 0) && (b[k][ii].p == i); ii--) b[k][ii].p = iii;
                }
                i = b[j][i].p;
        }
    #if BLAS
        cblas_daxpy(j, -mus, mu[j], 1, mu[k], 1);
    #else
        for (i = 0; i < j; i++) mu[k][i] -= mu[j][i]*mus;
    #endif

    }
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

/**
 * LLL variants
 */
void lll(lattice_t *lattice, int s, int z, DOUBLE quality, int deepinsert_blocksize) {
    DOUBLE **R;
    DOUBLE *beta;
    DOUBLE *N;
    DOUBLE **H;
    int r, bit_size;

    lllalloc(&R, &beta, &N, &H, s, z);
    bit_size = get_bit_size(lattice);
    r = lllHfp(lattice, R, beta, H, 0, 0, s, z, quality, deepinsert_blocksize, bit_size);
    lllfree(R, beta, N, H, s);

    return;
}

DOUBLE iteratedlll(lattice_t *lattice, int s, int z, int no_iterates, DOUBLE quality, int deepinsert_blocksize) {
    DOUBLE **R;
    DOUBLE *beta;
    DOUBLE *N;
    DOUBLE **H;
    int r, l, i, j, runs;
    int bit_size;
    coeff_t *swapvl;
    DOUBLE lD;

    lllalloc(&R, &beta, &N, &H, s, z);

    bit_size = get_bit_size(lattice);

    r = lllHfp(lattice, R, beta, H, 0, 0, s, z, quality, deepinsert_blocksize, bit_size);
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
            }
        }
        r = lllHfp(lattice, R, beta, H, 0, 0, s, z, quality, deepinsert_blocksize, bit_size);
        lD = log_potential(R, s, z);
        fprintf(stderr, "%d: log(D)= %f\n", runs, lD);
        fflush(stdout);
    }

    lllfree(R, beta, N, H, s);

    return lD;
}

DOUBLE block_reduce(lattice_t *lattice, int s, int z, int block_size, DOUBLE quality, int deepinsert_blocksize) {
    DOUBLE **R;
    DOUBLE *beta;
    DOUBLE *N;
    DOUBLE **H;
    DOUBLE lD;
    int start = 0, up, size, bit_size;
    coeff_t **basis_org;

    lllalloc(&R, &beta, &N, &H, s, z);
    bit_size = get_bit_size(lattice);

    //r = lllHfp(lattice, R, beta, H, 0, 0, s, z, quality, deepinsert_blocksize, bit_size);
    #if 1
    while (start < s) {
        fprintf(stderr, "Block reduce %d\n", start);
        up = start + block_size;
        up = (up > s) ? s : up;

        basis_org = lattice->basis;
        lattice->basis = &(lattice->basis[start]);
        size = (start + block_size > up) ? up - start : block_size;
        lllHfp(lattice, R, beta, H, 0, 0, size, z, quality, deepinsert_blocksize, bit_size);
        lattice->basis = basis_org;
        start += block_size;
    }
    #else
    while (start < s) {
        fprintf(stderr, "Block reduce %d\n", start);
        up = start + block_size;
        up = (up > s) ? s : up;

        lllHfp(lattice, R, beta, H, start, start, up, z, quality, deepinsert_blocksize, bit_size);
        start += block_size;
    }
    #endif
    //print_lattice(lattice);

    lD = log_potential(R, s, z);
    fprintf(stderr, "   log(D)= %f\n", lD);
    fflush(stderr);
    lllfree(R, beta, N, H, s);

    return lD;
}

void insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv) {
    coeff_t **b = lattice->basis;
    coeff_t *swapvl;
    int i, j, g;
    long q, ui;

    /* build new basis */
    for (j = 1; j <= z; j++)
        mpz_set_si(lattice->swap[j].c, 0);

    // Store new linear combination in lattice->swap
    for (i = start; i <= end; i++) {
        if (u[i] != 0) for(j = 1; j <= z; j++) {
            if (u[i] > 0) {
                mpz_addmul_ui(lattice->swap[j].c, b[i][j].c, u[i]);
            } else {
                mpz_submul_ui(lattice->swap[j].c, b[i][j].c, -u[i]);
            }
        }
    }
    g = end;
    while (u[g] == 0) g--;

    i = g - 1;
    while (labs(u[g]) > 1) {
        while (u[i] == 0) i--;
        q = (long)ROUND((1.0 * u[g]) / u[i]);
        ui = u[i];
        u[i] = u[g] - q*u[i];
        u[g] = ui;

        for (j = 1; j <= z; j++) {
            mpz_set(hv, b[g][j].c);
            mpz_mul_si(b[g][j].c, b[g][j].c, (long)q);
            mpz_add(b[g][j].c, b[g][j].c, b[i][j].c);
            mpz_set(b[i][j].c, hv);
        }
        coeffinit(b[g], z);
        coeffinit(b[i], z);
    }

    swapvl = b[g];
    for (i = g; i > start; i--)
        b[i] = b[i - 1];
    b[start] = lattice->swap;
    coeffinit(b[start], z);

    lattice->swap = swapvl;
    for (j = 1; j <= z; j++)
        mpz_set_si(lattice->swap[j].c, 0);
    coeffinit(lattice->swap, z);

}

/**
 * Blockwise Korkine Zolotareff reduction
 */
DOUBLE bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, int p) {
    DOUBLE **R, *h_beta, *N;
    DOUBLE **H;
    DOUBLE r_tt;
    DOUBLE new_cj;
    DOUBLE lD;

    static mpz_t hv;
    int zaehler;
    int h, i, last;
    int start_block, end_block;
    int bit_size = get_bit_size(lattice);

    long *u;

    mpz_init(hv);

    last = s - 1;    /* |last| points to the last nonzero vector of the lattice.*/
    if (last < 1) {
        printf("BKZ: the number of basis vectors is too small.\n");
        printf("Probably the number of rows is less or equal");
        printf(" to number of columns in the original system\n");
        printf("Maybe you have to increase c0 (the first parameter)!\n");

        mpz_clear(hv);
        return 0.0;
    }

    fprintf(stderr, "######### BKZ ########\n");
    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) {
        u[i] = 0;
    }

    lllalloc(&R, &h_beta, &N, &H, s, z);
    lllHfp(lattice, R, h_beta, H, 0, 0, s, z, delta, 10, bit_size); // delta

    start_block = zaehler = -1;
    //start_block = 0;
    while (zaehler < last) {
        start_block++;
        if (start_block == last)
            start_block = 0;

        end_block = start_block + beta - 1;
        end_block = (end_block < last) ? end_block : last;

        new_cj = enumerate(lattice, R, u, s, start_block, end_block, delta, p);
        //new_cj = sample(lattice, R, u, s, start_block, last);

        h = (end_block + 1 < last) ? end_block + 1 : last;

        r_tt = R[start_block][start_block];
        r_tt *= r_tt;
        if (delta * r_tt > new_cj) {
            fprintf(stderr, "enumerate successful %d %lf improvement: %lf\n",
                start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
            fflush(stderr);

            /* successful enumeration */
            insert_vector(lattice, u, start_block, end_block, z, hv);

            lllHfp(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, delta, 10, bit_size);
            //start_block = lllHfp(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, delta, 10, bit_size);
            //fprintf(stderr, "%d\n", start_block);

            zaehler = -1;
        } else {
            //fprintf(stderr, "enumerate: no improvement %d\n", zaehler);
            fflush(stderr);
            if (h > 0) {
                lllHfp(lattice, R, h_beta, H, h - 1, h - 1, h + 1, z, 0.0, -1, bit_size);
            }
            //start_block++;
            zaehler++;
        }
    } /* end of |while| */

    lD = log_potential(R, s-1, z);

    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stderr);
    lllfree(R, h_beta, N, H, s);
    free(u);
    mpz_clear(hv);

    return lD;
}

/**
 * Pruned Gauss-Enumeration.
 */
DOUBLE enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s,
                    int start_block, int end_block, DOUBLE improve_by, int p) {
    DOUBLE x;
    DOUBLE *y, *c;
    DOUBLE c_min;

    int i, j;
    int t, t_max;
    int found_improvement = 0;

    long *delta, *d, *v;
    DOUBLE *u_loc;
    int len, k;
    double alpha, radius;
    DOUBLE *lambda_min;
    int SCHNITT = 10;

    c = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    y = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    delta = (long*)calloc(s+1,sizeof(long));
    d = (long*)calloc(s+1,sizeof(long));
    v = (long*)calloc(s+1,sizeof(long));
    u_loc = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    lambda_min = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));

    len = end_block + 1 - start_block;
    for (i = start_block; i <= end_block + 1; i++) {
        c[i] = y[i] = 0.0;
        u_loc[i] = 0.0;
        v[i] = delta[i] = 0;
        d[i] = 1;
    }

    p = 0.25;
    if (end_block - start_block <= SCHNITT) {
        c_min = set_prune_const(R, start_block, end_block + 1, PRUNE_NO, 1.0);
    } else {
        //Hoerners version of the Gaussian volume heuristics.
        //hoerner(R, start_block, end_block + 1, p, eta);
        c_min = set_prune_const(R, start_block, end_block + 1, PRUNE_BKZ, p);
    }
    c_min *= improve_by;

    i = start_block;
    lambda_min[i] = R[i][i] * R[i][i];
    for (i = start_block + 1; i <= end_block; ++i) {
        x = R[i][i] * R[i][i];
        lambda_min[i] = (x < lambda_min[i-1]) ? x : lambda_min[i-1];
    }

    //t = t_max = end_block;
    for (t_max = start_block + 1; t_max <= end_block; t_max ++) {
        t = t_max;
        u_loc[t] = 1.0;

        while (t <= t_max) {
            handle_signals(lattice);

            x = (u_loc[t] + y[t]) * R[t][t];
            c[t] = c[t + 1] + x * x;

            if (len <= SCHNITT) {
                alpha = 1.0;
            } else {
                #if 0
                    alpha = 1.05 * (end_block + 1 - t) / len;
                #elif 1
                    k = (end_block + 1 - t);
                    if (k > 2 * len / 4) {
                        alpha = 1.0;
                    } else {
                        alpha = 0.6;
                        //alpha = k / len;
                    }
                #elif 0
                    alpha = 1.0;
                #else
                    p = 0.5;
                    k = (end_block + 1 - t);
                    if (k > len / 2) {
                        alpha = p * 2 * k / len;
                        alpha = 1.0;
                    } else {
                        alpha = 2 * p - 1 + 2 * k * (1 - p) / len;
                    }
                #endif
                alpha = (alpha < 1.0) ? alpha : 1.0;
            }
            //fprintf(stderr, "%d %d %d %lf\n", start_block, t, end_block, alpha);
            radius = alpha * c_min;

            #if 0
            if (t - start_block + 1 > 5) {
                x = lambda_min[t] * (t - start_block + 1) / 32;
                //fprintf(stderr, "> %lf\n", x);
                radius -= x;
            }
            #endif

            if (c[t] < radius - EPSILON) {
                if (t > start_block) {
                    // forward
                    t--;

                    #if BLAS
                        y[t] = cblas_ddot(t_max - t, &(u_loc[t+1]), 1, &(R[t+1][t]), lattice->num_rows);
                    #else
                        for (j = t + 1, y[t] = 0.0; j <= t_max; j++) {
                            y[t] += u_loc[j] * R[j][t];
                        }
                    #endif
                    y[t] /= R[t][t];

                    u_loc[t] = v[t] = (long)(ROUND(-y[t]));
                    delta[t] = 0;
                    d[t] = (v[t] > -y[t]) ? -1 : 1;

                    continue;
                } else {
                    // Found shorter vector
                    c_min = c[t];
                    for (i = start_block; i <= end_block; i++) {
                        u[i] = (long)round(u_loc[i]);
                        fprintf(stderr, "%ld ", u[i]);
                    }
                    fprintf(stderr, "\n");
                    found_improvement = 1;
                }
            } else {
                // back
                t++;
            }
            // next
            if (t < t_max) delta[t] *= -1.0;
            if (delta[t] * d[t] >= 0) delta[t] += d[t];
            u_loc[t] = v[t] + delta[t];
        }
    }

    free(c);
    free(y);
    free(delta);
    free(d);
    free(v);
    free(u_loc);
    free(lambda_min);

    if (!found_improvement) {
        c_min = R[start_block][start_block];
        c_min *= c_min;
    }
    return (c_min);
}

DOUBLE sample(lattice_t *lattice, DOUBLE **R, long *u, int s, int start_block, int end_block) {
    DOUBLE x;
    DOUBLE *y, *c;
    DOUBLE c_min;

    int i, j;
    int t, t_max;
    int found_improvement = 0;

    long *delta, *d, *v;
    DOUBLE *u_loc;
    double alpha, radius;

    c = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    y = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    delta = (long*)calloc(s+1,sizeof(long));
    d = (long*)calloc(s+1,sizeof(long));
    v = (long*)calloc(s+1,sizeof(long));
    u_loc = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));

    for (i = start_block; i <= end_block + 1; i++) {
        c[i] = y[i] = 0.0;
        u_loc[i] = 0.0;
        v[i] = delta[i] = 0;
        d[i] = 1;
    }

    radius = set_prune_const(R, start_block, end_block + 1, PRUNE_BKZ, 0.25);
    c_min = radius;

    //t = t_max = end_block;
    for (t_max = start_block + 1; t_max <= end_block; t_max ++) {
        t = t_max;
        u_loc[t] = 1.0;

        while (t <= t_max) {
            handle_signals(lattice);

            x = (u_loc[t] + y[t]) * R[t][t];
            c[t] = c[t + 1] + x * x;
            alpha = c_min;

            if (c[t] < alpha - EPSILON) {
                if (t > start_block) {
                    // forward
                    t--;

                    #if 0 //BLAS
                        y[t] = cblas_ddot(t_max - t, &(u_loc[t+1]), 1, &(R[t+1][t]), lattice->num_rows);
                    #else
                        for (j = t + 1, y[t] = 0.0; j <= t_max; j++) {
                            y[t] += u_loc[j] * R[j][t];
                        }
                    #endif
                    y[t] /= R[t][t];

                    u_loc[t] = v[t] = (long)(ROUND(-y[t]));
                    delta[t] = 0;
                    d[t] = (v[t] > -y[t]) ? -1 : 1;

                    continue;
                } else {
                    c_min = c[t];
                    for (i = start_block; i <= end_block; i++) {
                        u[i] = (long)round(u_loc[i]);
                        fprintf(stderr, "%ld ", u[i]);
                    }
                    fprintf(stderr, "\n");
                    found_improvement = 1;
                }
            } else {
                // back
                t++;
            }

            if (t < t_max - 10) {
                t++;
            } else if (delta[t] < 0) {
                t++;
            }

            // next
            if (t < t_max) delta[t] *= -1.0;
            if (delta[t] * d[t] >= 0) delta[t] += d[t];
            u_loc[t] = v[t] + delta[t];
        }
    }

    free(c);
    free(y);
    free(delta);
    free(d);
    free(v);
    free(u_loc);

    if (!found_improvement) {
        c_min = R[start_block][start_block];
        c_min *= c_min;
    }
    return (c_min);
}

/**
 * Exhaustive enumeration
*/
#define FINCKEPOHST 1

/**
 * Globals for enumeration
 */
#if FINCKEPOHST
    DOUBLE **muinv;
#endif
/*|mpz_t *upb,*lowb;|*/
long dual_bound_success;
DOUBLE dum1, dum2;

long only_zeros_no, only_zeros_success, hoelder_no, hoelder_success;
long hoelder2_success;
long cs_success;

DOUBLE explicit_enumeration(lattice_t *lattice, int columns, int rows) {
    /* local variables for |explicit_enumeration() */
    /*|__attribute((aligned(16)))|*/

    int level,level_max;
    int i,j,l;
    long loops;

    DOUBLE *y, *cs, *us;

    long *delta, *d, *eta;
    long *v;
    int *first_nonzero, *first_nonzero_in_column, *firstp;
    int bit_size;

    DOUBLE *N, **mu, *c, **w, **bd, **mu_trans;
    DOUBLE *bd_1norm;

    DOUBLE Fd, Fq, Fqeps;
    DOUBLE *dum;
    DOUBLE tmp;
    // coeff_t *swap_vec;

    //int isSideStep = 0;
    //DOUBLE stepWidth = 0.0;
    //DOUBLE olddum;

#if defined(FINCKEPOHST)
    DOUBLE *fipo;
    DOUBLE **dual_basis;
    DOUBLE *dual_bound;
#endif

    /* Vector to collect enumeration statistics */
    long nlow[1000];
    for (i=0;i<1000;i++) nlow[i] = 0;

    /* test the size of the basis */
    fprintf(stderr, "Dimension of solution space (k): %d compared to s-z+2: %d\n",
                columns, system_columns-system_rows+1+free_RHS);
    fflush(stderr);

    if (columns < system_columns - system_rows + 1 + free_RHS) {
        fprintf(stderr,"LLL didn't succeed in computing a basis of the kernel.\n");
        fprintf(stderr,"Please increase c0 (the first parameter)!\n");
        return 0;
    }

    /* allocate the memory for enumeration */
    lllalloc(&mu,&c,&N,&bd,columns,rows);

    us = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    cs = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    bd_1norm = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    y = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    delta = (long*)calloc(columns+1,sizeof(long));
    d = (long*)calloc(columns+1,sizeof(long));
    first_nonzero = (int*)calloc(rows,sizeof(int));
    first_nonzero_in_column = (int*)calloc(columns+rows+1,sizeof(int));
    if (first_nonzero_in_column == NULL)
        return(0);
    firstp = (int*)calloc(columns+1,sizeof(int));

    eta = (long*)calloc(columns+1,sizeof(long));
    v = (long*)calloc(columns+1,sizeof(long));
    w = (DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
    for (i = 0; i <= columns; i++)
        w[i] = (DOUBLE*)calloc(rows,sizeof(DOUBLE));

    mu_trans = (DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
    for (i = 0; i <= columns; i++)
        mu_trans[i]=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    dum = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));

#if FINCKEPOHST
    fipo = (DOUBLE*)calloc(columns+1, sizeof(DOUBLE));
    muinv = (DOUBLE**)calloc(columns, sizeof(DOUBLE*));
    for(i = 0; i < columns; ++i)
        muinv[i] = (DOUBLE*)calloc(rows, sizeof(DOUBLE));

    dual_basis = (DOUBLE**)calloc(columns+1, sizeof(DOUBLE*));
    for (i = 0; i <= columns; ++i) {
        dual_basis[i] = (DOUBLE*)calloc(rows, sizeof(DOUBLE));
    }
    dual_bound = (DOUBLE*)calloc(columns+1, sizeof(DOUBLE));

#endif

    /* initialize arrays */
    for (i = 0; i <= columns; i++) {
        cs[i] = y[i] = us[i] = 0.0;
        delta[i] = 0;
        v[i] = 0;
        eta[i] = d[i] = 1;
        for (l = 0; l < rows; l++)
            w[i][l] = 0.0;
    }

    bit_size = get_bit_size(lattice);

    /* count nonzero entries in the last rows(s) */
    if (free_RHS) {
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

#if 0
    @<sort lattice columns@>;
#endif

    /* set the simple pruning bounds */
    Fq = (DOUBLE)mpz_get_d(lattice->max_norm);
    Fd = (rows*Fq*Fq) * (1.0 + EPSILON);
    Fqeps = (1.0 + EPSILON) * Fq;        /* Used in prune() */
#if VERBOSE > 0
    fprintf(stderr, "Fq: %f\n", (double)Fq);
    fprintf(stderr, "Fd: %f\n", (double)Fd);
    fflush(stderr);
#endif

    /* orthogonalize the basis */
#if GIVENS
    givens(lattice, columns, rows, mu, bd, c);
#else
    gramschmidt(lattice, columns, rows, mu, bd, c);
#endif

    /* compute $mu^\top$, the transpose of $mu$. */
    for (i = 0; i < columns; i++)
        for (j = 0; j < columns; j++)
            mu_trans[j][i] = mu[i][j];

    /* Compute 1-norm of orthogonal basis */
    for (i = 0; i <= columns; ++i) {
        bd_1norm[i] = 0.0;
        for (j = 0; j < rows; ++j) {
            bd_1norm[i] += fabs(bd[i][j]);
        }
        bd_1norm[i] *= Fqeps / c[i];
    }

#if FINCKEPOHST
    /* determine Fincke-Pohst bounds */
    dual_bound_success = 0;
    inverse(mu, muinv, columns);

#if VERBOSE > -1
    fprintf(stderr, "Dual bounds:\n");
    fflush(stderr);
#endif

    /* Symmetric Fincke-Pohst */
    for (i = 0; i < columns; i++) {
        dum1 = 0.0;
        fipo[i] = 0.0;
        for (j = 0; j < rows; j++) {
            tmp = 0.0;
            for (l = i; l < columns; l++) {
                tmp += muinv[i][l] * bd[l][j] / c[l];
            }
            dual_basis[i][j] = tmp;
            fipo[i] += tmp * tmp;
            dum1 += fabs(tmp);
        }
        fipo[i] = SQRT(fipo[i] * Fd);
        dum1 =  fabs(dum1 * Fq) * (1.0 + EPSILON);
        if (dum1 < fipo[i]) {
            fipo[i] = dum1;
        }

#if VERBOSE > -1
        fprintf(stderr, "%0.3lf ", fipo[i]);
#endif
    }

    for (i = columns - 2; i >= 0; --i) {
        for (j = 0, tmp = 0.0; j < rows; j++) {
            dum1 = dual_basis[i][j] + dual_basis[i + 1][j];
            tmp += fabs(dum1);
        }
        dual_bound[i] = tmp * Fq * (1.0 + EPSILON);
    }

#if VERBOSE > -1
    fprintf(stderr, "\n\n");
    fflush(stderr);
#endif

#endif

    /* Remove trailing unnecessary columns. That means, columns
       whose corresponding Finke-Pohst bounds are equal to 0
       can be removed.
       This is important for the Selfdual Bent Functions Problems
     */
#if 1
    for (i = columns - 1; i >= 0; i--) {
        if (fipo[i]<0.9) {
            printf("DEL\n");
            columns--;
        } else {
            break;
        }
    }
#endif

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
    if (level<0) level = 0;
    level_max = level;
    us[level] = 1;
    v[level] = 1;
    only_zeros_no = only_zeros_success = 0;
    hoelder_no = hoelder_success = hoelder2_success = 0;
    cs_success = nosolutions = loops = 0;

    /* the loop of the exhaustive enumeration */
    do {
        /* increase loop counter */
        loops++;
        if ((lattice->LLL_params.stop_after_loops > 0) && (lattice->LLL_params.stop_after_loops <= loops))
            goto afterloop;

#if VERBOSE > -1
        if (loops % 100000000 ==0) {                 /*10000000*/
            printf("%ld loops, solutions: %ld ", loops, nosolutions);
#if FINCKEPOHST
            printf(", dual bounds: %ld ", dual_bound_success);
#endif
            printf("\n");
            fflush(stdout);
        }
#endif
        handle_signals(lattice);

        /* compute new |cs| */
        //olddum = dum[level];
        dum[level] = us[level] + y[level];
        cs[level] = cs[level+1] + dum[level]*dum[level]*c[level];

        if (cs[level] < Fd)  {
            /* Use (1, -1, 0, ...) as values in Hoelder pruning */
            if (fabs(dum[level]) > bd_1norm[level]) {
                ++hoelder2_success;
                goto step_back;
            }

#if FINCKEPOHST
            if (fabs(us[level]) > fipo[level]) {
                dual_bound_success++;
                goto side_step;
            }
#endif

            /*
            if (level < columns - 1 && fabs(us[level] + us[level + 1]) > dual_bound[level]) {
                fprintf(stderr, "%lf %lf\n", us[level] + us[columns -1], dual_bound[level]);
                fflush(stderr);
                goto side_step;
            }
            */

            /*
            if (isSideStep) {
                stepWidth = dum[level] - olddum;
                compute_w2(w, bd, stepWidth, level, rows);
            } else {
                compute_w(w, bd, dum[level], level, rows);
            }
            */
            compute_w(w, bd, dum[level], level, rows);

            if (level > 0) {
                /* not at a leave */
                i = prune_only_zeros(w, level, rows, Fq, first_nonzero_in_column, firstp,
                                     bd, y, us, columns);

                if (i < 0) {
                    goto step_back;
                } else if (i > 0) {
                    goto side_step;
                }

                ++hoelder_no;
                if (prune(w[level], cs[level], rows, Fqeps)) {
                    ++hoelder_success;
                    if (eta[level] == 1) {
                        goto step_back;
                    }
                    eta[level] = 1;
                    delta[level] *= -1;
                    if (delta[level]*d[level]>=0) delta[level] += d[level];
                    us[level] = v[level] + delta[level];
                } else {
                    level--;
                    delta[level] = eta[level] = 0;
                    y[level] = compute_y(mu_trans,us,level,level_max);
                    us[level] = v[level] = ROUND(-y[level]);
                    d[level] = (v[level]>-y[level]) ? -1 : 1;
                    //isSideStep = 0;
                }
            } else {
                /* at $|level|=0$ */
                if (final_test(w[0], rows, Fq, us, lattice, bit_size) == 1) {
                    print_solution(lattice, w[level], rows, Fq, us, columns);

                    if (lattice->LLL_params.stop_after_solutions > 0 &&
                        lattice->LLL_params.stop_after_solutions <= nosolutions)
                        goto afterloop;
                }
                goto side_step;


            }
        } else {
            cs_success++;
step_back:
            /* Up: we go to $|level|+1$. */
            nlow[level]++;
            level++;
            if (level_max<level) level_max = level;
side_step:
            /*
                Side step: the next value in the same level is
                chosen.
            */
            if (eta[level] == 0) {
                delta[level] *= -1;
                if (delta[level]*d[level] >= 0)
                    delta[level] += d[level];
            } else {
                delta[level] += (delta[level]*d[level]>=0) ? d[level] : -d[level] ;
            }
            us[level] = v[level] + delta[level];
            //isSideStep = 1;
        }
    } while (level<columns);
afterloop:

    /* final output */
    fprintf(stderr, "Prune_cs: %ld\n", cs_success);
    fprintf(stderr, "Prune_only_zeros: %ld of %ld\n", only_zeros_success, only_zeros_no);
    fprintf(stderr, "Prune_hoelder: %ld of %ld\n", hoelder_success, hoelder_no);
    fprintf(stderr, "Prune_hoelder interval: %ld\n", hoelder2_success);
#if FINCKEPOHST
    printf("Fincke-Pohst: %ld\n", dual_bound_success);
#endif
    printf("Loops: %ld\n",loops);

    if ((lattice->LLL_params.stop_after_solutions <= nosolutions && lattice->LLL_params.stop_after_solutions > 0) ||
            (lattice->LLL_params.stop_after_loops <= loops && lattice->LLL_params.stop_after_loops > 0 )) {
        printf("Stopped after number of solutions: %ld\n", nosolutions);

        if (lattice->LLL_params.silent)
            print_num_solutions(nosolutions);
        if ((lattice->LLL_params.stop_after_loops <= loops && lattice->LLL_params.stop_after_loops > 0)) {
            exit(10);
        } else {
            exit(9);
        }
    } else {
        printf("Total number of solutions: %ld\n", nosolutions);
    }
    printf("\n");
    fflush(stdout);

    /* free allocated memory for enumeration */
    free(us);
    free(cs);
    free(bd_1norm);
    free(y);
    free(delta);
    free(d);
    free(first_nonzero);
    free(first_nonzero_in_column);
    free(firstp);

    free(eta);
    free(v);
    for (l = 0; l <= columns; l++) free(w[l]);
    free(w);
    free(original_columns);

#if FINCKEPOHST
    free(fipo);
    for (l = 0; l < columns; l++) free(muinv[l]);
    free(muinv);

    for (l = 0; l <= columns; l++) free(dual_basis[l]);
    free(dual_basis);
    free(dual_bound);
#endif

    lllfree(mu, c, N, bd, columns);
    for (l = 0; l < columns; l++) free(mu_trans[l]);
    free(mu_trans);

    return 1;
}

DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max) {
    #if BLAS
        return cblas_ddot(level_max-level, &(us[level+1]), 1, &(mu_trans[level][level+1]), 1);
    #else
        int i;
        DOUBLE dum;
        i = level_max;
        dum = 0.0;
        while (i >= level + 1) {
            dum += mu_trans[level][i]*us[i];
            i--;
        }
        return dum;
    #endif
}

void compute_w2(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows) {
    #if BLAS
        cblas_daxpy(rows, alpha, bd[level], 1, w[level], 1);
    #else
        int i;
        for (i = rows - 1; i >= 0; --i) {
            /*
            DOUBLE x, y;
            x = w[level][i] + alpha * bd[level][i];
            y = w[level+1][i] + beta * bd[level][i];

            if (fabs(x - y) > 0.0000000001) {
                fprintf(stderr,"ERROR w2 level=%d: %lf %lf %lf %lf \n", level, x, y, alpha, beta);
                fprintf(stderr,"                  :%lf %lf %lf %lf\n", old, ys, u, bd[level][i]);
                fprintf(stderr,"                  :%lf\n", (x - y)/ bd[level][i]);
                fflush(stderr);
            }


            w[level][i] = x;
            */
            w[level][i] += alpha * bd[level][i];
        }
    #endif

    return;
}

void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows) {
    #if BLAS
        cblas_dcopy(rows, w[level+1], 1, w[level], 1);
        cblas_daxpy(rows, alpha, bd[level], 1, w[level], 1);
    #else
        int i;

        i = rows - 1;
        while (i >= 0) {
            w[level][i] = w[level+1][i] + alpha * bd[level][i];
            i--;
        }
    #endif

    return;
}

void gramschmidt(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c) {
    int i, l, j;
    DOUBLE dum;

    for (i = 0; i < columns; i++) {
        for (l = 0; l < rows; l++) bd[i][l] = (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l));
        for (j = 0; j < i; j++) {
            dum = 0.0;
            for (l = 0; l < rows; l++) dum += (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l)) * bd[j][l];
            mu[i][j] = dum / c[j];
            for (l = 0; l < rows; l++) bd[i][l] -= mu[i][j]*bd[j][l];
        }

        c[i] = scalarproductfp(bd[i], bd[i], rows);
        #if VERBOSE > 0
            printf("%lf ",(double)c[i]);
        #endif
    }
    #if VERBOSE > 0
        printf("\n\n");
        fflush(stdout);
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
            printf("%6.3f ",(double)c[i]);
            if (i>0 && i%15==0) printf("\n");
        #endif
    }
    #if VERBOSE > -1
        printf("\n\n");
        fflush(stdout);
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

    // If the involved numbers are to big,
    // an exact test is done.
    if (bit_size < 27) {
        return 1;
    }
    for (i = 0; i < rows; i++) {
        if (!iszeroone) {
            if (mpz_cmp_si(upperbounds[i], 0) != 0) {
                mpz_divexact(soltest_upfac, upperbounds_max, upperbounds[i]);
            } else {
                mpz_set(soltest_upfac, upperbounds_max);
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
        mpz_divexact(soltest_u, soltest_u, max_norm_initial);
        mpz_divexact(soltest_u, soltest_u, soltest_upfac);
        mpz_divexact_ui(soltest_u, soltest_u, denom);
        mpz_abs(soltest_u, soltest_u);
        if (!iszeroone && (mpz_cmp_si(soltest_u, 0) < 0 || mpz_cmp(soltest_u, upperbounds[i]) > 0) ) {
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

int prune_only_zeros(DOUBLE **w, int level, int rows, DOUBLE Fq,
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

        if (iszeroone) {
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
    upper = rows-1-free_RHS;
    if (free_RHS && fabs(w[upper])>Fq*(1+EPSILON)) {
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
        end = (lattice->cut_after == -1) ? no_original_columns : lattice->cut_after;
        for (j = 0; j < end; j++) {
            if (original_columns[j] == 0) {
                mpz_set_si(soltest_u,0);
            } else {
                if (!iszeroone) {
                    if (mpz_cmp_si(upperbounds[i],0)!=0) {
                        mpz_divexact(soltest_upfac,upperbounds_max,upperbounds[i]);
                    } else {
                        mpz_set(soltest_upfac,upperbounds_max);
                    }
                }
                mpz_set_si(soltest_u,0);
                for (k=0;k<columns;k++) {
                    if (ROUND(us[k])>0) {
                        mpz_addmul_ui(soltest_u,get_entry(lattice->basis, k, i), ROUND(us[k]));
                    } else {
                        mpz_submul_ui(soltest_u,get_entry(lattice->basis, k, i), -ROUND(us[k]));
                    }
                }
                mpz_sub(soltest_u, soltest_u, soltest_s);
                mpz_divexact(soltest_u, soltest_u, max_norm_initial);
                mpz_divexact(soltest_u, soltest_u, soltest_upfac);
                mpz_divexact_ui(soltest_u, soltest_u, denom);
                mpz_abs(soltest_u, soltest_u);

                i++;
            }
            mpz_out_str(NULL, 10, soltest_u);
            fflush(stdout);
            mpz_out_str(fp, 10, soltest_u);

            /* Meanwhile, all solution vectors are written with separating blanks. */
            /*|if (!iszeroone) { }|*/
            printf(" ");
            fprintf(fp, " ");
        }
        if (free_RHS) {
            mpz_set_d(soltest_u,ROUND(w[i]));
            mpz_divexact(soltest_u,soltest_u,max_up);
            mpz_abs(soltest_u,soltest_u);
            printf(" L = ");
            mpz_out_str(NULL,10,soltest_u);
        }
        printf("\n");
        fflush(stdout);
        fprintf(fp, "\n");
        fflush(fp);
    }

    nosolutions++;
    if (nosolutions%10000==0) {
        printf("%ld\n",nosolutions);
        fflush(stdout);
    }

    return 1;
}

void stop_program_sig(int sig) {
    if (sig != SIGALRM)
       return;

    printf("Stopped after SIGALRM, number of solutions: %ld\n", nosolutions);
    if (!SILENT)
        print_num_solutions(nosolutions);

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

void handle_signals(lattice_t *lattice) {
    if (PRINT_REQUIRED) {
        print_lattice(lattice);
        PRINT_REQUIRED = 0;
    }
    if (DUMP_REQUIRED) {
        dump_lattice(lattice);
        DUMP_REQUIRED = 0;
    }
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

/*
    An estimate on gamma_1(L[low, up]), excluding up
    lambda_1(L) = (det(L) / V_n(1))^(1/n)

    V_n(1) = pi^(n/2) / Gamma(n/2 + 1)
    Gamma(n/2 + 1) = sqrt(pi n) (n/2)^(n/2)*e^(-n/2)
 */
DOUBLE GH(DOUBLE **R, int low, int up) {
    int i, n, k;
    DOUBLE log_det, V1;
    static DOUBLE pi = 3.141592653589793238462643383;
    //static DOUBLE e = 2.718281828459045235360287471352662497757247093;

    for (i = low, log_det = 0.0; i < up; ++i) {
        log_det += log(R[i][i] * R[i][i]);
    }
    log_det *= 0.5;
    n = up - low;

    // Exact formulae for unit ball volume
    if (n % 2 == 0) {
        k = n / 2;
        for (i = 1, V1 = 1.0; i <= k; i++) {
            V1 *= pi / i;
        }
    } else {
        k = (n - 1)/ 2;
        for (i = 0, V1 = 1.0 / pi; i <= k; i++) {
            V1 *= 2.0 * pi / (2*i + 1);
        }
    }
    V1 = exp(log(V1) / n);

    return exp(log_det / n) / V1;
}

/*
    Hoerners version of the Gaussian volume heuristics.
*/
void hoerner(DOUBLE **R, int low, int up, double p, DOUBLE *eta) {
    int i;
    static DOUBLE pi = 3.141592653589793238462643383;
    static DOUBLE e = 2.718281828459045235360287471352662497757247093;
    DOUBLE c, x;
    int t_up;

    c = R[low][low];
    x = log(c * c);
    for (i = low + 1; i < up; i++) {
        t_up = i - low;
        eta[i] = 0.5 * t_up * exp((log(pi * t_up) - 2.0 * p * log(2.0) + x) / t_up ) / (pi * e);
        //fprintf(stderr, "::: %lf \n", eta[i]);
        if (i < up - 1) {
            c = R[i][i];
            x += log(c * c);
        }
    }
}

DOUBLE set_prune_const(DOUBLE **R, int low, int up, int prune_type, DOUBLE p) {
    DOUBLE gh, gh1;
    DOUBLE alpha;

    alpha = 1.05;

    gh1 = R[low][low];
    gh1 *= gh1;

    if (prune_type == PRUNE_BKZ) {
        gh = GH(R, low, up);
        gh *= gh;
        gh *= alpha;
    } else {
        gh = gh1; // prune_type == PRUNE_NO
    }

    fprintf(stderr, ">>>> %d %lf %lf\n", up - low, gh1, gh);
    fflush(stderr);

    return (gh <= gh1) ? gh : gh1;
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
    mpz_out_str(NULL, 10, upperbounds_max);
    printf("\n\n[");
    for (i = 0; i < lattice->num_rows - 1; i++) {
        mpz_out_str(NULL, 10, upperbounds[i]);
        printf(" ");
    }
    printf("]\n");
    fflush(stdout);

    return;
}
