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
long num_solutions;
int nboundvars;
int level_max;

mpz_t soltest_u;
mpz_t soltest_s;
mpz_t soltest_upfac;

static FILE* fp;

long diophant(lgs_t *LGS, lattice_t *lattice, FILE* solfile, int restart, char *restart_filename) {

    int i, j;
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
    lattice->basis_long = (long**)calloc(lattice->num_cols + ADDITIONAL_COLS, sizeof(long*));
    for (j = 0; j < lattice->num_cols/* + ADDITIONAL_COLS*/; j++) {
        lattice->basis[j] = (coeff_t*)calloc(lattice->num_rows + 1, sizeof(coeff_t));
        lattice->basis_long[j] = (long*)calloc(lattice->num_rows, sizeof(long));
        for (i = 0; i <= lattice->num_rows; i++)
            mpz_init(lattice->basis[j][i].c);
    }
    lattice->swap = (coeff_t*)calloc(lattice->num_rows + 1, sizeof(coeff_t));
    lattice->swap_long = (long*)calloc(lattice->num_rows + 1, sizeof(long));
    for (i = 0; i <= lattice->num_rows; i++) {
        mpz_init(lattice->swap[i].c);
    }
    lattice->work_on_long = FALSE;

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
        lll(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_LOW, CLASSIC_LLL);

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

    /**
     * third reduction
     */
    fprintf(stderr, "\n"); fflush(stderr);

    if (lattice->LLL_params.iterate) {
        iteratedlll(lattice, lattice->num_cols, lattice->num_rows, lattice->LLL_params.iterate_no, LLLCONST_HIGH, POT_LLL);
    } else {
        //shufflelattice(lattice);

        #if 0
        i = 0;
        do {
            lD = lDnew;

            //shufflelattice(lattice);
            #if 0
            if (i % 2 == 0) {
                lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGHER,
                        lattice->LLL_params.bkz.beta, lattice->LLL_params.bkz.p,
                        solutiontest, solutiontest_long);
                fprintf(stderr, "BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);

                dump_lattice(lattice);
            } else {
                lDnew = dual_bkz(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGHER,
                            lattice->LLL_params.bkz.beta, lattice->LLL_params.bkz.p,
                            solutiontest);
                fprintf(stderr, "Dual BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);
            }
            fflush(stderr);
            #endif

            lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGHER,
                        lattice->LLL_params.bkz.beta, lattice->LLL_params.bkz.p,
                        solutiontest, solutiontest_long);
            fprintf(stderr, "BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);

            i++;
        } while (i < 1 && fabs(lDnew - lD) > 0.01);
        #else
            lDnew = self_dual_bkz(lattice, lattice->num_cols, lattice->num_rows, LLLCONST_HIGHER,
                    lattice->LLL_params.bkz.beta, lattice->LLL_params.bkz.p,
                    solutiontest);
        #endif
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
        mpz_out_str(stderr,10,soltest_u);
        fprintf(stderr, " ");
        if (lattice->LLL_params.stop_after_solutions == 1) {
            mpz_out_str(fp,10,soltest_u);
            fprintf(fp," ");
        }
    }
    if (free_RHS) {
        mpz_divexact(soltest_u, get_entry(lattice->basis, position, up), max_up);
        mpz_divexact(soltest_u, soltest_u, lattice->LLL_params.scalelastlinefactor);
        mpz_abs(soltest_u,soltest_u);
        fprintf(stderr, " L = ");
        mpz_out_str(stderr,10,soltest_u);
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
        mpz_out_str(stderr,10,soltest_u);
        fprintf(stderr, " ");
        if (lattice->LLL_params.stop_after_solutions == 1) {
            mpz_out_str(fp,10,soltest_u);
            fprintf(fp," ");
        }
    }
    if (free_RHS) {
        mpz_divexact(soltest_u, get_entry(lattice->basis, position, up), max_up);
        mpz_divexact(soltest_u, soltest_u, lattice->LLL_params.scalelastlinefactor);
        mpz_abs(soltest_u,soltest_u);
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
    DOUBLE **R;
    DOUBLE *beta;
    DOUBLE *N;
    DOUBLE **H;
    int r, bit_size;

    lllalloc(&R, &beta, &N, &H, s, z);
    bit_size = get_bit_size(lattice);
    r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, solutiontest);
    lllfree(R, beta, N, H, s);

    return;
}

DOUBLE iteratedlll(lattice_t *lattice, int s, int z, int no_iterates, DOUBLE quality, int reduction_type) {
    DOUBLE **R;
    DOUBLE *beta;
    DOUBLE *N;
    DOUBLE **H;
    int r, l, i, j, runs;
    int bit_size;
    coeff_t *swapvl;
    long *swap;
    DOUBLE lD;

    lllalloc(&R, &beta, &N, &H, s, z);

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
    lllfree(R, beta, N, H, s);

    return lD;
}

DOUBLE block_reduce(lattice_t *lattice, int s, int z, int block_size, DOUBLE quality, int reduction_type) {
    DOUBLE **R;
    DOUBLE *beta;
    DOUBLE *N;
    DOUBLE **H;
    DOUBLE lD;
    int start = 0, up, size, bit_size;
    coeff_t **basis_org;

    lllalloc(&R, &beta, &N, &H, s, z);
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
    lllfree(R, beta, N, H, s);

    return lD;
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

typedef struct {
    DOUBLE diff;
    DOUBLE cs;
    DOUBLE l1;
    DOUBLE y;
    int num;
    int pos;
    DOUBLE us;
    DOUBLE* w;
} enum_node_t;

typedef struct {
    int pos;
    int num;
    enum_node_t* nodes;
} enum_level_t;

int cmpfunc (const void* a, const void* b) {
    DOUBLE ea = ((enum_node_t*)a)->diff;
    DOUBLE eb = ((enum_node_t*)b)->diff;
    //fprintf(stderr, "%lf %lf\n", ((enum_node_t*)a)->diff, ((enum_node_t*)b)->diff);

    return (ea - eb);
}

void mysort(enum_node_t *data, int len) {
    enum_node_t swap;
    int i, j, pos;

    for (i = 0 ; i < len -1; i++) {
        pos = i;

        for (j = i + 1 ; j < len ; j++) {
            if (data[pos].diff > data[j].diff) {
                pos = j;
            }
        }
        if (pos != i) {
            swap = data[i];
            data[i] = data[pos];
            data[pos] = swap;
        }
    }
}

typedef struct {
    long loops;

    long *delta;
    long *d;
    long *eta;
    long *v;

    DOUBLE *y;
    DOUBLE *cs;
    DOUBLE *us;

    DOUBLE *dum;
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
    zigzag->dum = (DOUBLE*)calloc(c, sizeof(DOUBLE));

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

int enumLevel(enum_level_t* enum_data, zigzag_t* zigzag, lattice_t* lattice,
                DOUBLE** bd, DOUBLE* c,
                DOUBLE Fd, DOUBLE Fqeps, DOUBLE Fq,
                DOUBLE* bd_1norm, DOUBLE* fipo,
                int level, int rows, int columns, int bit_size) {

    DOUBLE olddum, stepWidth;
    int i, j, isSideStep;

    enum_data[level].num = 0;
    isSideStep = FALSE;
    do {
        /* increase loop counter */
        zigzag->loops++;
        if ((lattice->LLL_params.stop_after_loops > 0) &&
            (lattice->LLL_params.stop_after_loops <= zigzag->loops)) {
            //goto afterloop;
            return -1;
        }

        #if VERBOSE > -1
        if (zigzag->loops % 100000000 ==0) {                 /*10000000*/
            fprintf(stderr, "%ld loops, solutions: %ld",
                zigzag->loops, num_solutions);
            #if FINCKEPOHST
                fprintf(stderr, ", dual bounds: %ld ", dual_bound_success);
            #endif
            fprintf(stderr, "\n");
            fflush(stderr);
        }
        #endif

        handle_signals(lattice, NULL);

        /* compute new |cs| */
        olddum = zigzag->dum[level];
        zigzag->dum[level] = zigzag->us[level] + zigzag->y[level];
        zigzag->cs[level] = zigzag->cs[level+1] + zigzag->dum[level] * zigzag->dum[level] * c[level];

        if (zigzag->cs[level] < Fd)  {
            /* Use (1, -1, 0, ...) as values in Hoelder pruning */
            if (fabs(zigzag->dum[level]) > bd_1norm[level]) {
                //++hoelder2_success;
                // goto recurse;
                break;
            }

            #if FINCKEPOHST
            if (fabs(zigzag->us[level]) > fipo[level]) {
                //dual_bound_success++;
                isSideStep = FALSE;
                goto side_step;
            }
            #endif

            // if (isSideStep) {
            //     stepWidth = zigzag->dum[level] - olddum;
            //     compute_w2(zigzag->w, bd, stepWidth, level, rows);
            // } else {
                compute_w(zigzag->w, bd, zigzag->dum[level], level, rows);
            // }

            if (level >= 0) {
                /* not at a leave */
                // i = prune_only_zeros(w, level, rows, Fq, first_nonzero_in_column, firstp,
                //                      bd, y, us, columns);
                // if (i < 0) {
                //    //goto step_back;
                //    goto recurse;
                // } else if (i > 0) {
                //     goto side_step;
                // }

                //++hoelder_no;
                if (prune(zigzag->w[level], zigzag->cs[level], rows, Fqeps)) {
                    //++hoelder_success;
                    if (zigzag->eta[level] == 1) {
                        //goto step_back;
                        //goto recurse;
                        break;
                    }
                    zigzag->eta[level] = 1;
                    zigzag->delta[level] *= -1;
                    if (zigzag->delta[level] * zigzag->d[level] >= 0) {
                        zigzag->delta[level] += zigzag->d[level];
                    }
                    zigzag->us[level] = zigzag->v[level] + zigzag->delta[level];
                    isSideStep = TRUE;
                } else {
                    i = enum_data[level].num;
                    enum_data[level].nodes[i].us = zigzag->us[level];
                    enum_data[level].nodes[i].cs = zigzag->cs[level];
                    for (j = 0; j < rows; j++) {
                        enum_data[level].nodes[i].w[j] = zigzag->w[level][j];
                    }
                    enum_data[level].nodes[i].l1 = cblas_dasum(rows, zigzag->w[level], 1);
                    enum_data[level].nodes[i].y = zigzag->y[level];
                    enum_data[level].nodes[i].diff =
                        //-(Fqeps * enum_data[level].nodes[i].l1 - zigzag->cs[level]);
                        //cblas_dasum(rows, zigzag->w[level], 1);
                        zigzag->cs[level];
                        //zigzag->cs[level] / cblas_dasum(rows, zigzag->w[level], 1);
                        //-(Fqeps *cblas_dasum(rows, zigzag->w[level], 1) - zigzag->cs[level])
                        //    / zigzag->cs[level];
                    enum_data[level].num++;

                    isSideStep = TRUE;
                    goto side_step;
                }
            }
        } else {
            //cs_success++;
            return 0;

side_step:
            /*
                Side step: the next value in the same level is
                chosen.
            */
            if (zigzag->eta[level] == 0) {
                zigzag->delta[level] *= -1;
                if (zigzag->delta[level] * zigzag->d[level] >= 0) {
                    zigzag->delta[level] += zigzag->d[level];
                }
            } else {
                zigzag->delta[level] += zigzag->d[level] *
                    ((zigzag->delta[level] * zigzag->d[level] >= 0) ? 1: -1);
            }
            zigzag->us[level] = zigzag->v[level] + zigzag->delta[level];
        }

    } while (TRUE);

    return 0;
}

int dfs(enum_level_t* enum_data, zigzag_t* zigzag, lattice_t* lattice,
                DOUBLE** bd, DOUBLE* c,
                DOUBLE Fd, DOUBLE Fqeps, DOUBLE Fq,
                DOUBLE* bd_1norm, DOUBLE* fipo,
                int level, int rows, int columns, int bit_size, DOUBLE** mu_trans) {

    int j, i;
    DOUBLE s;

    if (-1 == enumLevel(enum_data, zigzag, lattice,
            bd, c, Fd, Fqeps, Fq, bd_1norm, fipo,
            level, rows, columns, bit_size)) {

        return -1;
    }

    //qsort(enum_data[level].nodes, enum_data[level].num, sizeof(enum_node_t), cmpfunc);
    //mysort(enum_data[level].nodes, enum_data[level].num);

    #if 0
        if (enum_data[level].num > 0) {
            fprintf(stderr, "level %d %0.2lf %0.2lf\n",
                level, enum_data[level].nodes[0].us, -zigzag->y[level]);
        }
    #endif

    for (enum_data[level].pos = 0;
            enum_data[level].pos < enum_data[level].num;
            enum_data[level].pos++) {

        // if ((level == 78 || level == 71 || level == 67) && enum_data[level].pos == 0) {
        //      continue;
        // }
        zigzag->us[level] = enum_data[level].nodes[enum_data[level].pos].us;
        zigzag->cs[level] = enum_data[level].nodes[enum_data[level].pos].cs;

        // zigzag->w[level] = enum_data[level].nodes[enum_data[level].pos].w;
        for (j = 0; j < rows; j++) {
            zigzag->w[level][j] = enum_data[level].nodes[enum_data[level].pos].w[j];
        }

        if (level == 0) {
            // Solution found
            if (final_test(zigzag->w[0], rows, Fq, zigzag->us, lattice, bit_size) == 1) {
                print_solution(lattice, zigzag->w[level], rows, Fq, zigzag->us, columns);

                for (j = columns - 1 ; j >= 0; j--) {
                    //if (1 || zigzag->us[j] != ROUND(-zigzag->y[j])) {
                    if (enum_data[j].pos > 0) {
                        fprintf(stderr, "================== ");
                    }
                        fprintf(stderr, "%d: %0.lf %0.lf %d of %d:\n",
                            j, zigzag->us[j], ROUND(-zigzag->y[j]),
                            enum_data[j].pos, enum_data[j].num
                        );
                        for (i = 0;
                             i <= enum_data[j].num - 1;
                             i++) {
                            s = enum_data[j].nodes[i].y + enum_data[j].nodes[i].us;
                            fprintf(stderr, "\t%.0lf\t%lf\t%lf\t%lf\t dum=%lf\n",
                                enum_data[j].nodes[i].us,
                                enum_data[j].nodes[i].cs,
                                enum_data[j].nodes[i].l1,
                                s * s * c[j],
                                s
                            );
                        }
                }

                if (lattice->LLL_params.stop_after_solutions > 0 &&
                    lattice->LLL_params.stop_after_solutions <= num_solutions)
                    return -1; //goto afterloop;
            }
            continue;
        }

        level--;

        zigzag->delta[level] = zigzag->eta[level] = 0;
        enum_data[level].pos = 0;
        zigzag->y[level] = compute_y(mu_trans, zigzag->us, level, level_max);
        zigzag->us[level] = zigzag->v[level] = ROUND(-zigzag->y[level]);
        zigzag->d[level] = (zigzag->v[level] > -zigzag->y[level]) ? -1 : 1;

        if (-1 == dfs(enum_data, zigzag, lattice,
                bd, c, Fd, Fqeps, Fq, bd_1norm, fipo,
                level, rows, columns, bit_size, mu_trans)) {
            return -1;
        }

        level++;
    }

    level++;
    if (level >= columns) {
        // We are done, let's leave the loop.
        //break;
        return 1;
    } else if (level > level_max) {
        level_max = level;
    }
    return 0;
}

int lds(enum_level_t* enum_data, zigzag_t* zigzag, lattice_t* lattice,
                DOUBLE** bd, DOUBLE* c,
                DOUBLE Fd, DOUBLE Fqeps, DOUBLE Fq,
                DOUBLE* bd_1norm, DOUBLE* fipo,
                int level, int rows, int columns, int bit_size, DOUBLE** mu_trans,
                int lds_k, int lds_l) {

    int j, i;
    int result;
    int start, end, pos, do_left_branch, p;
    int next_lds_k;
    DOUBLE s;

    if (-1 == enumLevel(enum_data, zigzag, lattice,
            bd, c, Fd, Fqeps, Fq, bd_1norm, fipo,
            level, rows, columns, bit_size)) {

        return -1;
    }

    #if 0
    if (enum_data[level].num == 0) {
        fprintf(stderr, "-----------\n"); fflush(stderr);
    }
    #endif

    // if (lds_k == 0) {
    //     start = 0;
    //     if (level < lds_l) {
    //         end = enum_data[level].num;
    //     } else {
    //         end = (1 <= enum_data[level].num) ? 1 : enum_data[level].num;
    //     }
    // } else if (lds_k == 1) {
    //     start = 0;
    //     end = (2 <= enum_data[level].num) ? 2 : enum_data[level].num; // min
    // } else {
    //     start = 0;
    //     end = (2 <= enum_data[level].num) ? 2 : enum_data[level].num; // min
    // }

    start = 1;
    do_left_branch = 1;
    if (level <= lds_l) {
        end = enum_data[level].num;
    } else {
        if (level - lds_l <= lds_k) {   // depth <= k -> no left branch
            start = 1;
            do_left_branch = 0;
        }
        if (lds_k > 0) {
            end = (lds_k < enum_data[level].num) ? lds_k + 1 : enum_data[level].num;
        } else {
            // left-branches only
            end = 1;
        }
    }

    #if 0
    if (level >= lds_l) {
        if (level == level_max) {
            fprintf(stderr, "=============================================\n");
        }
        fprintf(stderr, "LDS %d %d from %d to %d\n", level, lds_k, start, end); fflush(stderr);
    }
    #endif
    #if 0
    if (level < 25 && lds_k > 0) {
        fprintf(stderr, "*");
        fflush(stderr);
    }
    #endif

    // for (pos = start; pos < end && pos < enum_data[level].num; pos++) {
    for (pos = start; pos <= enum_data[level].num; pos++) {
        if (pos >= end &&
            !(do_left_branch && pos == enum_data[level].num)) {
                continue;
            }
        p = pos % enum_data[level].num;
        enum_data[level].pos = p;

        zigzag->us[level] = enum_data[level].nodes[p].us;
        zigzag->cs[level] = enum_data[level].nodes[p].cs;
        // zigzag->w[level] = enum_data[level].nodes[p].w;
        for (j = 0; j < rows; j++) {
            zigzag->w[level][j] = enum_data[level].nodes[p].w[j];
        }

        if (level == 0) {
            // Solution found
            if (final_test(zigzag->w[0], rows, Fq, zigzag->us, lattice, bit_size) == 1) {
                print_solution(lattice, zigzag->w[level], rows, Fq, zigzag->us, columns);

                for (j = columns - 1 ; j >= 0; j--) {
                    //if (1 || zigzag->us[j] != ROUND(-zigzag->y[j])) {
                    if (p > 0) {
                        fprintf(stderr, "================== ");
                    }
                        fprintf(stderr, "%d: %d of %d\n",
                            j, //zigzag->us[j], ROUND(-zigzag->y[j]),
                            enum_data[j].pos, enum_data[j].num //, fabs( enum_data[j].nodes[0].us + zigzag->y[j])
                        );
                        for (i = 0; 0 && i <= enum_data[j].num - 1; i++) {
                            s = enum_data[j].nodes[i].y + enum_data[j].nodes[i].us;
                            fprintf(stderr, "\t%.0lf\t%lf\t%lf\t%lf\t dum=%lf\n",
                                enum_data[j].nodes[i].us,
                                enum_data[j].nodes[i].cs,
                                enum_data[j].nodes[i].l1,
                                s * s * c[j],
                                s
                            );
                        }
                }

                if (lattice->LLL_params.stop_after_solutions > 0 &&
                    lattice->LLL_params.stop_after_solutions <= num_solutions)
                    return -1; //goto afterloop;
            }
            continue;
        }

        level--;

        zigzag->delta[level] = zigzag->eta[level] = 0;
        enum_data[level].pos = 0;
        zigzag->y[level] = compute_y(mu_trans, zigzag->us, level, level_max);
        zigzag->us[level] = zigzag->v[level] = ROUND(-zigzag->y[level]);
        zigzag->d[level] = (zigzag->v[level] > -zigzag->y[level]) ? -1 : 1;

        next_lds_k = lds_k;
        if (level > lds_l) {
            // we are in ILDS mode
            if (p == 0) {
                // depth > k, left branch
                next_lds_k = lds_k;
            } else if (lds_k > 0) {
                next_lds_k = (lds_k > p) ? lds_k - p : 0;
            }
        }

        result = lds(enum_data, zigzag, lattice,
                bd, c, Fd, Fqeps, Fq, bd_1norm, fipo,
                level, rows, columns, bit_size, mu_trans, next_lds_k, lds_l);

        if (result == -1) {
            return result;
        }

        level++;
    }

    level++;
    if (level >= columns) {
        // We are done, let's leave the loop.
        //break;
        return 1;
    } else if (level > level_max) {
        level_max = level;
    }
    return 0;
}

DOUBLE explicit_enumeration(lattice_t *lattice, int columns, int rows) {
    /* local variables for |explicit_enumeration() */
    /*|__attribute((aligned(16)))|*/

    int level;
    //level_max;
    int i, j, l, k;
    int result;

    // DOUBLE *y, *cs, *us;
    // long *delta, *d, *eta;
    // long *v;
    // DOUBLE *dum;
    zigzag_t zigzag;
    enum_level_t* enum_data;

    DOUBLE *bd_1norm;
    int *first_nonzero, *first_nonzero_in_column, *firstp;
    int bit_size;

    DOUBLE *N, **mu, *c, **bd, **mu_trans;

    DOUBLE Fd, Fq, Fqeps;
    DOUBLE tmp;
    // coeff_t *swap_vec;

    int isSideStep = FALSE;
    DOUBLE stepWidth = 0.0;
    DOUBLE olddum;

    #if defined(FINCKEPOHST)
    DOUBLE *fipo;
    DOUBLE **dual_basis;
    DOUBLE *dual_bound;
    #endif


    /* Vector to collect enumeration statistics */
    long nlow[1000];
    for (i = 0; i < 1000; i++) nlow[i] = 0;

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
    bd_1norm = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));

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
    // dum = (DOUBLE*)calloc(columns+1, sizeof(DOUBLE));

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
            fprintf(stderr, "DEL\n");
            columns--;
        } else {
            break;
        }
    }
#endif

    /* New strategy */
    enum_data = (enum_level_t*)calloc(columns+1, sizeof(enum_level_t));
    for (i = 0; i <= columns; i++) {
        enum_data[i].nodes = (enum_node_t*)calloc(20 * ((int)(fipo[i]) + 1), sizeof(enum_node_t));
        for (j = 0; j < 20 * ((int)(fipo[i]) + 1); j++) {
            enum_data[i].nodes[j].w = (DOUBLE*)calloc(rows, sizeof(DOUBLE));
        }
        enum_data[i].num = 0;
        enum_data[i].pos = 0;
    }

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

    level = columns - 1;
    level_max = level;

    allocateZigzag(&zigzag, columns, rows, level);

    only_zeros_no = only_zeros_success = 0;
    hoelder_no = hoelder_success = hoelder2_success = 0;
    cs_success = 0;

    /* the loop of the exhaustive enumeration */
    #if 1
        //for (i = 0; i <= columns / 2; i++) {
        for (k = 0; k <= columns; k++) {
            fprintf(stderr, "lds_k=%d\n", k); fflush(stderr);
            result = lds(enum_data, &zigzag, lattice,
                bd, c, Fd, Fqeps, Fq, bd_1norm, fipo,
                level, rows, columns, bit_size, mu_trans,
                k, 0);
                //k, columns - columns / 2);
            //if (k == 1) break;
            if (result  == -1) {
                fprintf(stderr, "solution for lds_k=%d\n\n", k);
                break;
            }
        }
    #else
        dfs(enum_data, &zigzag, lattice,
            bd, c, Fd, Fqeps, Fq, bd_1norm, fipo,
            level, rows, columns, bit_size, mu_trans);
    #endif

//afterloop:
    /* final output */
    fprintf(stderr, "Prune_cs: %ld\n", cs_success);
    fprintf(stderr, "Prune_only_zeros: %ld of %ld\n", only_zeros_success, only_zeros_no);
    fprintf(stderr, "Prune_hoelder: %ld of %ld\n", hoelder_success, hoelder_no);
    fprintf(stderr, "Prune_hoelder interval: %ld\n", hoelder2_success);
#if FINCKEPOHST
    fprintf(stderr, "Dual bounds: %ld\n", dual_bound_success);
#endif
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

#if FINCKEPOHST
    // free(fipo);
    // for (l = 0; l < columns; l++) free(muinv[l]);
    // free(muinv);
    //
    // for (l = 0; l <= columns; l++) free(dual_basis[l]);
    // free(dual_basis);
    // free(dual_bound);
#endif

    // lllfree(mu, c, N, bd, columns);
    // for (l = 0; l < columns; l++) free(mu_trans[l]);
    // free(mu_trans);

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
