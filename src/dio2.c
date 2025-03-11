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
#include "enum.h"

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

/**
 * global variables
 */
// Defined in datastruct.h
extern long num_solutions;

// Defined in sd2.c
extern bool do_dump;

// Defined in enum.h
extern solution_t solution;

//mpz_t max_norm_initial;
//mpz_t max_up;
mpz_t dummy;
mpz_t lastlines_factor;

long diophant(lgs_t *LGS, lattice_t *lattice, FILE* solfile, int restart, char *restart_filename) {

    int i;
    int block_size;
    double lD = 0.0;
    double lDnew = 0.0;
    mpz_t *swap_vec;

    /**
     * Initialize some globals
     */
    mpz_init(lastlines_factor);
    mpz_init(solution.u);
    mpz_init(solution.s);
    mpz_init_set_ui(solution.upfac, 1);

    #if BLAS
        fprintf(stderr, "Use OpenBLAS: %s\n", openblas_get_config());
        openblas_set_num_threads(4);
    #endif

    if (!preprocess(LGS)) {
        fprintf(stderr, "Preprocess found contradiction\n");
        fprintf(stderr, "Total number of solutions: 0\n\n");
        exit(EXIT_NOT_SOLVABLE);
    }
    fprintf(stderr, "Num variables after preprocess: %d\n", LGS->num_cols);

    /*
     * Generate lattice from LGS and do scaling.
     * Allocate memory for lattice struct.
     */
    lgs_to_lattice(LGS, lattice);
    #if 0
        fprintf(stderr, "After scaling\n");
        print_lattice(lattice, stderr);
    #endif

    /**
     * open solution file
     */
    solution.fp = solfile;
    if (lattice->LLL_params.silent) fprintf(solution.fp, "SILENT\n");
    fflush(solution.fp);


    #if 1 // Do reduction

        /**
         * Rotate last lattice column to the first position
         */
        swap_vec = lattice->basis[lattice->num_cols-1];
        for (i = lattice->num_cols - 1; i > 0; i--) {
            lattice->basis[i] = lattice->basis[i - 1];
        }
        lattice->basis[0] = swap_vec;
        //shufflelattice(lattice);

        #if IS_USED
           fprintf(stderr, "After permute\n");
           print_lattice(lattice, stderr);
        #endif

        /*
         * ------ First reduction -----
         */
        mpz_set_ui(lastlines_factor, 1);
        fprintf(stderr, "\n"); fflush(stderr);
        if (!restart) {
            int org_cols = lattice->num_cols;  // Kernel vectors are swapped to the end.
            lll(lattice, lattice->num_cols, lattice->num_rows, lattice->LLL_params.lll.delta_low, DEEP_LLL/*KERNEL_LLL*//* DEEP_LLL */);
            lattice->num_cols = org_cols;

            #if IS_USED
                fprintf(stderr, "After first reduction\n");
                print_lattice(lattice, stderr);
                fprintf(stderr, "max norm ");
                mpz_out_str(stderr, 10, lattice->max_norm);
                fprintf(stderr, "\n");
                // exit(0);
            #endif

            /*
             * Cut the lattice
             */
            fprintf(stderr, "Before cutting\n");
            if (cutlattice(lattice)) {
                fprintf(stderr, "First reduction successful\n"); fflush(stderr);
            } else {
                fprintf(stderr, "First reduction not successful\n"); fflush(stderr);
                return 0;
            }
            #if IS_USED
                fprintf(stderr, "After cutting\n");
                print_lattice(lattice, stderr);
            #endif
            // exit(-1);

            #if 1
                /**
                 * ------ Second reduction -----
                 */
                // shufflelattice(lattice);
                mpz_set_ui(lastlines_factor, 1);
                lll(lattice, lattice->num_cols, lattice->num_rows, lattice->LLL_params.lll.delta_med, POT_LLL);
                lll(lattice, lattice->num_cols, lattice->num_rows, lattice->LLL_params.lll.delta_high, POT_LLL);
                fprintf(stderr, "Second reduction successful\n"); fflush(stderr);
            #endif


        } else {
            load_lattice(lattice, restart_filename);
            //print_lattice(lattice, stderr);
        }

        #if IS_USED
            fprintf(stderr, "After second reduction\n");
            print_lattice(lattice, stderr);
        #endif

        #if TRUE
        /**
         * ------ Third reduction ------
         */

        // Scale last rows
        mpz_set(lastlines_factor, lattice->LLL_params.scalelastlinefactor);
        for (i = 0; i < lattice->num_cols; i++) {
            mpz_mul(lattice->basis[i][lattice->num_rows - 1], lattice->basis[i][lattice->num_rows - 1], lastlines_factor);
            lattice->basis_long[i][lattice->num_rows - 1] *= mpz_get_si(lastlines_factor);
        }
        if (lattice->free_RHS) {
            for (i = 0; i < lattice->num_cols; i++) {
                mpz_mul(lattice->basis[i][lattice->num_rows - 2], lattice->basis[i][lattice->num_rows - 2], lastlines_factor);
                lattice->basis_long[i][lattice->num_rows - 2] *= mpz_get_si(lastlines_factor);
            }
        }

        //
        // Iterate or bkz
        //
        if (lattice->LLL_params.type == ITERATE) {
            // ---------- Iterated LLL
            iteratedlll(lattice, lattice->num_cols, lattice->num_rows,
                    lattice->LLL_params.iterate_no,
                    lattice->LLL_params.lll.delta_high, POT_LLL);
        } else if (lattice->LLL_params.type == PROGBKZ) {
            // ---------- Progressive BKZ
            for (block_size = 4; block_size <= lattice->LLL_params.bkz.beta; block_size += 4) {
                lD = lDnew;
                shufflelattice(lattice);
                lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows,
                            lattice->LLL_params.lll.delta_higher,
                            block_size, lattice->LLL_params.bkz.p,
                            ENUM_BLOCK, 30,
                            solutiontest, solutiontest_long);
                fprintf(stderr, "BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);
            }
        } else {
            // ---------- BKZ
            // shufflelattice(lattice);
            lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows,
                        lattice->LLL_params.lll.delta_higher,
                        lattice->LLL_params.bkz.beta, lattice->LLL_params.bkz.p,
                        ENUM_BLOCK, 10,
                        solutiontest, solutiontest_long);
            fprintf(stderr, "BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);

            #if IS_USED
                lD = lDnew;
                shufflelattice(lattice);
                lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows,
                            lattice->LLL_params.lll.delta_higher,
                            10000, lattice->LLL_params.bkz.p,
                            ENUM_LDS_FULL, 20,
                            solutiontest, solutiontest_long);
                fprintf(stderr, "LDS improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);
            #endif

            #if IS_USED
                for (i = 1; i < 10; i++) {
                    lD = lDnew;
                    shufflelattice(lattice);
                    lDnew = bkz(lattice, lattice->num_cols, lattice->num_rows,
                                lattice->LLL_params.lll.delta_higher,
                                lattice->LLL_params.bkz.beta + i, lattice->LLL_params.bkz.p,
                                ENUM_BLOCK, 200,
                                solutiontest, solutiontest_long);
                    fprintf(stderr, "BKZ improvement: %0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);
                }
            #endif
        }

        #if GSA_OUT==TRUE
            print_gsa(lattice->decomp.R, lattice->num_cols, 0);
        #endif
        fprintf(stderr, "Third reduction successful\n"); fflush(stderr);

        if (do_dump) {
            dump_lattice(lattice);
        }

        /*
         * Undo scaling of last rows
         */
        for (i = 0; i < lattice->num_cols; i++) {
            mpz_divexact(lattice->basis[i][lattice->num_rows - 1], lattice->basis[i][lattice->num_rows - 1], lastlines_factor);
            lattice->basis_long[i][lattice->num_rows - 1] /= mpz_get_si(lastlines_factor);
        }
        if (lattice->free_RHS) {
            for (i = 0; i < lattice->num_cols; i++) {
                mpz_divexact(lattice->basis[i][lattice->num_rows - 2], lattice->basis[i][lattice->num_rows - 2], lastlines_factor);
                lattice->basis_long[i][lattice->num_rows - 2] /= mpz_get_si(lastlines_factor);
            }
        }

        // print_lattice(lattice, stderr);

        #endif // End of third reduction
    #else
        read_NTL_lattice();
    #endif // Do reduction

    if (lattice->LLL_params.print_ntl) {
        fprintf(stderr, "Print lattice for NTL and exit\n");
        print_NTL_lattice(lattice);   /* Version for the NTL output */
        return 0;
    }

    if (lattice->LLL_params.kernel) {
        fprintf(stderr, "Print kernel of the system to stdout\n");
        print_kernel(lattice);
        return 0;
    }


    /**
     * explicit enumeration
     */
    fprintf(stderr, "\n"); fflush(stderr);
    num_solutions = exhaustive_enumeration(lattice);

    if (lattice->LLL_params.silent)
        print_num_solutions(num_solutions);

    /**
     * Free lattice memory
     */
    free_lattice(lattice);

    // for (i = 0; i < lattice->num_rows; i++) {
    //     mpz_clear(lattice->swap[i]);
    // }
    // free(lattice->swap);
    // mpz_clear(lattice->matrix_factor);
    // mpz_clear(lattice->max_norm);

    mpz_clear(lastlines_factor);
    mpz_clear(solution.u);
    mpz_clear(solution.s);
    mpz_clear(solution.upfac);

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
int cutlattice(lattice_t *lattice) {
    int j, i, flag;
    int all_zero;
    size_t rows_aligned;

    /**
     * Delete unnecessary columns
     */
    j = 0;
    do {
        all_zero = 1;
        for (i = 0; i < lattice->lgs_rows; ++i) {
            if (mpz_sgn(lattice->basis[j][i]) != 0) {
                all_zero = 0;
                break;
            }
        }
        if (all_zero) {
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
        if (mpz_sgn(get_entry(lattice->basis, i, lattice->num_rows - 1)) != 0) {
            flag = 1;
            break;
        }

    if (flag == 0) {
        fprintf(stderr, "Non-homogenous solution not possible.\n"); fflush(stderr);
        exit(EXIT_NOT_SOLVABLE);

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
    rows_aligned = DO_ALIGN(lattice->num_rows * sizeof(double)) / sizeof(double);
    for (i = 1; i < lattice->num_cols; i++) {
        lattice->decomp.mu[i] = (double*)(lattice->decomp.mu[0] + i * rows_aligned);
    }
    int m = (lattice->num_rows > lattice->num_cols) ? lattice->num_rows : lattice->num_cols;
    for (i = 1; i < m; i++) {
        lattice->decomp.bd[i] = (double*)(lattice->decomp.bd[0] + i * rows_aligned);
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
    //mpz_divexact(solution.s, get_entry(lattice->basis, position,
    //    lattice->num_rows - 1), lattice->LLL_params.scalelastlinefactor);
    mpz_set(solution.s, get_entry(lattice->basis, position, lattice->num_rows - 1));

    /* write a solution with blanks */
    i = low;

    end = (lattice->cut_after == -1) ? lattice->num_cols_org : lattice->cut_after;

    for (j = 0; j < end; j++) {
        if (lattice->original_cols[j] == 0) {
            mpz_set_si(solution.u, 0);
        } else {
            if (!lattice->is_zero_one) {
                if (mpz_cmp_si(lattice->upperbounds[i-low], 0) != 0) {
                    mpz_divexact(solution.upfac, lattice->upperbounds_max, lattice->upperbounds[i - low]);
                } else {
                    mpz_set(solution.upfac, lattice->upperbounds_max);
                }
            }
            mpz_set(solution.u, get_entry(lattice->basis, position, i));
            mpz_sub(solution.u, solution.u, solution.s);
            mpz_divexact(solution.u, solution.u, lattice->max_norm_initial);
            mpz_divexact(solution.u, solution.u, solution.upfac);
            mpz_divexact_ui(solution.u, solution.u, lattice->denom);
            mpz_abs(solution.u, solution.u);
            i++;
        }
        mpz_out_str(stderr, 10, solution.u);
        fprintf(stderr, " ");
        if (lattice->LLL_params.stop_after_solutions == 1) {
            mpz_out_str(NULL, 10, solution.u);
            mpz_out_str(solution.fp, 10, solution.u);
            printf(" ");
            fprintf(solution.fp," ");
        }
    }
    if (lattice->free_RHS) {
        mpz_divexact(solution.u, get_entry(lattice->basis, position, up), lattice->max_up);
        // mpz_divexact(solution.u, solution.u, lattice->LLL_params.scalelastlinefactor);
        mpz_abs(solution.u, solution.u);
        fprintf(stderr, " L = ");
        mpz_out_str(NULL, 10, solution.u);
        mpz_out_str(stderr, 10, solution.u);
    }
    fprintf(stderr, " !!\n");
    fflush(stderr);

    /* test if one solution is enough */
    if (lattice->LLL_params.stop_after_solutions == 1) {
        printf("\n");
        fprintf(solution.fp,"\n");

        fprintf(stderr, "Stopped in phase 1 after finding a random solution\n");
        exit(EXIT_RANDOM_SOLUTION);
    }

    return 1;
}

int solutiontest_long(lattice_t *lattice, int position) {
    int i, j;
    int low, up;
    int end;

    #if IS_USED
    int is_good = 1;
    for (j = 0; j < lattice->num_rows; ++j) {
        if (labs(lattice->basis_long[position][j]) != 1) {
            is_good = 0;
            break;
        }
    }
    if (is_good) {
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SOLUTION\n");
    }
    #endif

    for (j = 0; j < lattice->num_rows; ++j) {
        mpz_set_si(lattice->basis[position][j], lattice->basis_long[position][j]);
    }
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
    // mpz_divexact(solution.s, get_entry(lattice->basis, position, lattice->num_rows-1), lattice->LLL_params.scalelastlinefactor);
    mpz_set(solution.s, get_entry(lattice->basis, position, lattice->num_rows - 1));

    /* write a solution with blanks */
    i = low;
    end = (lattice->cut_after == -1) ? lattice->num_cols_org : lattice->cut_after;
    for (j = 0; j < end; j++) {
        if (lattice->original_cols[j] == 0) {
            mpz_set_si(solution.u,0);
        } else {
            if (lattice->is_zero_one) {
                if (mpz_cmp_si(lattice->upperbounds[i - low], 0) != 0) {
                    mpz_divexact(solution.upfac, lattice->upperbounds_max, lattice->upperbounds[i-low]);
                } else {
                    mpz_set(solution.upfac, lattice->upperbounds_max);
                }
            }
            mpz_set(solution.u,get_entry(lattice->basis, position, i));
            mpz_sub(solution.u, solution.u, solution.s);
            mpz_divexact(solution.u, solution.u, lattice->max_norm_initial);
            mpz_divexact(solution.u, solution.u, solution.upfac);
            mpz_divexact_ui(solution.u, solution.u, lattice->denom);
            mpz_abs(solution.u, solution.u);
            i++;
        }
        mpz_out_str(stderr,10, solution.u);
        fprintf(stderr, " ");
        if (lattice->LLL_params.stop_after_solutions == 1) {
            mpz_out_str(NULL, 10, solution.u);
            mpz_out_str(solution.fp, 10, solution.u);
            printf(" ");
            fprintf(solution.fp, " ");
        }
    }
    if (lattice->free_RHS) {
        mpz_divexact(solution.u, get_entry(lattice->basis, position, up), lattice->max_up);
        // mpz_divexact(solution.u, solution.u, lattice->LLL_params.scalelastlinefactor);
        mpz_abs(solution.u, solution.u);
        fprintf(stderr, " L = ");
        mpz_out_str(stderr, 10, solution.u);
    }
    fprintf(stderr, " ||\n");
    fflush(stderr);

    /* Test if one solution is enough */
    if (lattice->LLL_params.stop_after_solutions == 1) {
        printf("\n");
        fprintf(solution.fp,"\n");

        fprintf(stderr, "Stopped in phase 1 after finding a random solution\n");
        exit(EXIT_RANDOM_SOLUTION);
    }

    return 1;
}

/**
 * LLL variants
 */
void lll(lattice_t *lattice, int s, int z, double quality, int reduction_type) {
    double **R = lattice->decomp.R;
    double *beta = lattice->decomp.c;
    //double *N = lattice->decomp.N;
    double **H = lattice->decomp.H;
    // int r;
    int bit_size;

    bit_size = get_bit_size(lattice);
    // Ignore return value of lllH
    lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, WORDLEN_MPZ, solutiontest);

    return;
}

double iteratedlll(lattice_t *lattice, int s, int z, int no_iterates, double quality, int reduction_type) {
    double **R = lattice->decomp.R;
    double *beta = lattice->decomp.c;
    //double *N = lattice->decomp.N;
    double **H = lattice->decomp.H;
    int r, i, j, runs;
    int bit_size;
    mpz_t *swapvl;
    long *swap;
    double lD;

    bit_size = get_bit_size(lattice);

    if (bit_size < 32) {
        copy_lattice_to_long(lattice);
        // r = lllH_long(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, solutiontest_long);
        r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, WORDLEN_LONG, solutiontest);
    } else {
        r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, WORDLEN_MPZ, solutiontest);
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
            r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, WORDLEN_LONG, solutiontest_long);
        } else {
            r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, WORDLEN_MPZ, solutiontest);
        }
        lD = log_potential(R, s, z);
        fprintf(stderr, "%d: log(D)= %f\n", runs, lD);
        fflush(stderr);
        #if GSA_OUT==TRUE
            print_gsa(R, s, 0);
        #endif

    }
    if (bit_size < 32) {
        copy_lattice_to_mpz(lattice);
    }
    return lD;
}

double block_reduce(lattice_t *lattice, int s, int z, int block_size, double quality, int reduction_type) {
    double **R = lattice->decomp.R;
    double *beta = lattice->decomp.c;
    double **H = lattice->decomp.H;

    double lD;
    int start = 0, up, size, bit_size;
    mpz_t **basis_org;

    bit_size = get_bit_size(lattice);

    //r = lllH(lattice, R, beta, H, 0, 0, s, z, quality, reduction_type, bit_size, WORDLEN_MPZ);
    while (start < s) {
        fprintf(stderr, "Block reduce %d\n", start);
        up = start + block_size;
        up = (up > s) ? s : up;
        basis_org = lattice->basis;
        lattice->basis = &(lattice->basis[start]);
        size = (start + block_size > up) ? up - start : block_size;
        lllH(lattice, R, beta, H, 0, 0, size, z, quality, reduction_type, bit_size, WORDLEN_MPZ, solutiontest);
        lattice->basis = basis_org;
        start += block_size;
    }
    //print_lattice(lattice, stderr);

    lD = log_potential(R, s, z);
    fprintf(stderr, "   log(D)= %f\n", lD);
    fflush(stderr);

    return lD;
}

void stop_program_sig(int sig) {
    if (sig != SIGALRM)
       return;

    fprintf(stderr, "Stopped after SIGALRM, number of solutions: %ld\n", num_solutions);
    if (!SILENT)
        print_num_solutions(num_solutions);

    exit(EXIT_MAX_SIGALRM);
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
    int i, j;

    fprintf(stderr, "%d x %d\n", lattice->num_cols, lattice->num_rows);
    // printf("%d %d\n", lattice->num_cols, lattice->num_rows);
    printf("\n[");
    for (i = 0; i < lattice->num_cols; i++) {
        printf("[");
        for (j = 0; j < lattice->num_rows; j++) {
            mpz_out_str(NULL, 10, lattice->basis[i][j]);
            printf(" ");
        }
        printf("]");
        printf("\n");
    }
    printf("]\n");
    fflush(stdout);

    if (!lattice->is_zero_one) {
        printf("\n");

        //printf("%d ", lattice->num_cols - 1);
        mpz_out_str(NULL, 10, lattice->upperbounds_max);
        printf("\n\n[");
        for (i = 0; i < lattice->num_rows - 1; i++) {
            fprintf(stderr, "%d\n", i);
            mpz_out_str(NULL, 10, lattice->upperbounds[i]);
            printf(" ");
        }
        printf("]\n");
        fflush(stdout);
    }

    return;
}

/**
 * Output the kernel of the system of equations to stdout.
 * It is required that there are no upper bounds, i.e.
 * the appended matrix has 2s in the diagonal.
 */
void print_kernel(lattice_t *lattice) {
    int i, j,
        rank = 0;
    mpz_t q;

    for (i = 0; i < lattice->num_cols; i++) {
        if (mpz_cmp_si(lattice->basis[i][lattice->num_rows - 1], 0) == 0) {
            rank++;
        }
    }

    if (rank != lattice->num_cols - 1) {
        fprintf(stderr, "Kernel computation failed, please increase parameter -scalelastline\n");
        exit(1);
    }

    fprintf(stderr, "%d x %d\n", rank, lattice->num_rows - 1);
    printf("%d %d %d\n", rank, lattice->num_rows - 1, 1);

    mpz_init(q);
    for (i = 0; i < lattice->num_cols; i++) {
        if (mpz_cmp_si(lattice->basis[i][lattice->num_rows - 1], 0) != 0) {
            continue;
        }
        for (j = 0; j < lattice->num_rows - 1; j++) {
            mpz_divexact_ui(q, lattice->basis[i][j], 2);
            printf("%ld ", mpz_get_si(q));
        }
        printf(" 0\n");
    }
    fflush(stdout);
    mpz_clear(q);

    return;
}
