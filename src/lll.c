#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "const.h"
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
#else
    #define BLAS 0
#endif

/**
 * @brief LLL reduction algorithms of a basis consisting of
 * basis vectors of type *mpz (GMP integers).
 * It is possible to choose
 * between the reduction types POT_LLL, CLASSIC_LLL or deep insert.
 *
 * @param {lattice_t} lattice
 * @param {DOUBLE**} R Matrix R of the QR decomposition of the matrix consisting of the
 *        lattice vectors from position [start...up[ (both inclusive)
 * @param {DOUBLE*} beta
 * @param {DOUBLE**} H Matrix containing Householder columns
 * @param {int} start Position of the basis vector in the array of lattice vectors
 *              at which the reduction will start. If less than "low", "low" will be taken.
 * @param {int} low Position of the first basis vector in the array of lattice vectors.
 *              LLL will not go below this number.
 * @param {int} up First position of basis vector in the array of lattice vectors
 *              which will not be visited.
 * @param {int} z Number of "rows", i.e. length of the basis vectors.
 * @param {DOUBLE} delta Reduction quality (0.5 <= delta <= 1)
 * @param {int} reduction_type reduction type: CLASSIC_LLL, POT_LLL, else: deep insert
 * @param {int} bit_size
 * @param {function*} solutiontest
 * @return int The lowest position of lattice vectors which the algorithm has visited.
 */
int lllH(lattice_t *lattice, DOUBLE **R, DOUBLE *beta, DOUBLE **H,
            int start, int low, int up, int z,
            DOUBLE delta, int reduction_type,
            int bit_size, int word_len,
            int (*solutiontest)(lattice_t *lattice, int k)) {

    mpz_t **b = lattice->basis;
    long **bl = lattice->basis_long;

    int i, j, k;
    DOUBLE norm;
    DOUBLE theta, eta;

    DOUBLE mus;
    mpz_t musvl;
    mpz_t *swapvl;
    long *swap;

    DOUBLE r_new, r_act;
    DOUBLE pot, pot_max;
    int pot_idx;
    int insert_pos, lowest_pos, k_max;
    int deep_size;
    int redo_tricol = 0;
    int max_tricols = 0;
    int stop_tricol = 10;
    int count_tricols = 0;

    // fprintf(stderr, "delta=%lf\n", delta);
    #if VERBOSE > 1
        int counter = 0;
    #endif
    fprintf(stderr, "-------------------------- Do LLLH %d -------------------------------------------------------------\n", word_len);

    lattice->work_on_long = FALSE;

    eta = ETACONST;
    if (bit_size > 100) {
        eta = 0.52;
        theta = 0.50;
    } else if (bit_size > 55) {
        eta = 0.51;
        theta = 0.2;
    } else if (bit_size > 30) {
        eta = 0.505;
        theta = 0.1;
    } else {
        eta = ETACONST;
        theta = 0.0;
        // eta = 0.52;
        // theta = 0.3;
    }

    #if VERBOSE > 1
        fprintf(stderr, "LLLH: ");
        if (reduction_type == POT_LLL) {
            fprintf(stderr, "use PotLLL");
        } else if (reduction_type == CLASSIC_LLL) {
            fprintf(stderr, "use classic LLL");
        } else {
            fprintf(stderr, "use deepinsert %d", reduction_type);
        }
        fprintf(stderr, ", delta=%0.3lf", delta);
        fprintf(stderr, ", eta=%0.3lf", eta);
        fprintf(stderr, ", max bits: %d\n", bit_size);
    #endif

    if (word_len == WORDLEN_MPZ) {
        mpz_init(musvl);
    }

    /* Test for trivial cases. */
    if ((z <= 1) || (up <= 1)) {
        fprintf(stderr, "Wrong dimensions in lllH\n");
        fflush(stderr);
        return(0);
    }

    k = (start >= low) ? start : low;
    k_max = k;
    lowest_pos = k;

    /* The main loop */
    while (k < up) {
        if (k < lowest_pos) {
            lowest_pos = k;
        }
        #if VERBOSE > 1
            if (counter > 0 && (counter % 10000) == 0) {
                fprintf(stderr, "LLL: %d k:%d\n", counter, k);
                fflush(stderr);
            }
            counter++;
        #endif
        handle_signals(lattice, R);

        // #if FALSE
        // // Look ahead
        // i = householder_column(b, R, H, beta, k, s, z, bit_size);
        // if (i > k) {
        //     swapvl = b[i];
        //     for (j = i; j > k; --j) {
        //         b[j] = b[j - 1];
        //     }
        //     b[k] = swapvl;

        //     //fprintf(stderr, "GET %d from %d\n", k, i);
        // }
        // #endif

        count_tricols = 0;
    again:
        if (FALSE && count_tricols > 0) {
            fprintf(stderr, "\nBefore householder\n");
            if (word_len == WORDLEN_MPZ) {
                check_precision(b[k], R[k], z, k);
            }
        }

        // Apply Householder vectors to column k of R  and
        // determine the Householder vector for column k
        if (word_len == WORDLEN_MPZ) {
            i = householder_column(b, R, H, beta, k, k + 1, z, bit_size);
        } else {
            i = householder_column_long(bl, R, H, beta, k, k + 1, z, bit_size);
        }

        // if (fabs(R[k][k]) < 1.0e-12) {
        //     goto swap_zero_vector;
        // }
        if (FALSE && k > 0) {
            fprintf(stderr, "\nBefore: R[%d]: ", k);
            for (j = 0; j <= k; j++) {
                fprintf(stderr, " %0.20lf", R[k][j]);
            }
            fprintf(stderr, "\n");
            fprintf(stderr, "mu: %0.20lf\n", R[k][k-1] / R[k-1][k-1]);
        }

        if (FALSE && count_tricols > 0) {
            fprintf(stderr, "After householder\n");
            if (word_len == WORDLEN_MPZ) {
                check_precision(b[k], R[k], z, k);
            }
        }

        if (FALSE) fprintf(stderr, "k=%d\n", k);
        
        redo_tricol = 0;
        /* Size reduction of $b_k$ */
        for (j = k - 1; j >= low; j--) {
            /**
             * Subtract suitable multiple of $b_j$ from $b_k$.
             *
             * Size reduction according to Schnorr.
             * Lazy size reduction, see Stehle, "Floating-point LLL: theoretical and practical aspects"
             */
            if (fabs(R[k][j]) > eta * fabs(R[j][j]) + theta * fabs(R[k][k])) {
                mus = ROUND(R[k][j] / R[j][j]);
                redo_tricol = 1;
                // fprintf(stderr, "%0.2lf ", R[k][j] / R[j][j]);

                /* set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$ */
                if (word_len == WORDLEN_MPZ) {
                    mpz_set_d(musvl, mus);
                    size_reduction(b, R, musvl, mus, k, j, z);
                } else {
                    size_reduction_long(bl, R, (long)mus, mus, k, j, z);
                }
                (*solutiontest)(lattice, k);
            }
        }
        // fprintf(stderr, "\n");

        if (FALSE && redo_tricol > 0) {
            fprintf(stderr, "\n");
            fprintf(stderr, "redo=%d, tricol=%d \n", redo_tricol, count_tricols+1);
            if (word_len == WORDLEN_MPZ) {
                // check_precision(b[k - 1], R[k - 1], z, k - 1);
                check_precision(b[k], R[k], z, k);
            }
        }

        if (redo_tricol) {
            count_tricols++;
            max_tricols = (count_tricols > max_tricols) ? count_tricols : max_tricols;
            if (count_tricols > stop_tricol) {
                fprintf(stderr, "LLLH: too much tricol iterations (%d), k=%d, bit size: %d\n", stop_tricol, k, bit_size);

                // for (i = 0; i < k; i++) {
                //     fprintf(stderr, "%0.0lf ", R[i][i]);
                // }
                // fprintf(stderr, "\n");

                exit(EXIT_ERR_NUMERIC);
            }
            goto again;
        }

        if (FALSE && k > 0) {
            fprintf(stderr, "After: R[%d]: ", k);
            for (j = 0; j <= k; j++) {
                fprintf(stderr, " %0.20lf", R[k][j]);
            }
            fprintf(stderr, "\n");
            fprintf(stderr, "mu: %0.20lf\n", R[k][k-1] / R[k-1][k-1]);
        }

        /*
            Before going to step 4 we test if $b_k$ is linear dependent.
            If we find a linear dependent vector $b_k$,
            we shift b_k to the last column of the
            matrix and restart lllH with s = s-1.
        */
        // #if BLAS
        //     norm = cblas_dnrm2(k + 1, R[k], 1);
        // #else
        //     for (j = 0, norm = 0.0; j <= k; ++j) {
        //         norm += R[k][j] * R[k][j];
        //     }
        //     norm = SQRT(norm);
        // #endif
        // TODO
        // It is enough to check R[k][k]
        norm = hiprec_norm_l2(R[k], k + 1);

        if (norm != norm || norm < 0.5) {  // nan or < 0.5
    swap_zero_vector:
            // print_lattice(lattice);
            // fprintf(stderr, "lllH: Zero vector at %d\n", k);
            // fflush(stderr);
            // TODO
            // Check where the 0-vector is positioned to
            // lattice->num_cols - 1 or lattice->num_cols
            if (word_len == WORDLEN_MPZ) {
                swapvl = b[k];
                for (i = k + 1; i < lattice->num_cols; i++) {
                    b[i-1] = b[i];
                }
                b[lattice->num_cols - 1] = swapvl;
            } else {
                swap = bl[k];
                for (i = k + 1; i < lattice->num_cols; i++) {
                    bl[i-1] = bl[i];
                }
                bl[lattice->num_cols - 1] = swap;
            }

            up--;
            lattice->num_cols--;

            k = 0;
            continue;
        }

        // If delta == 0, only size reduction is done
        if (delta == 0.0) {
            k++;
            continue;
        }

        /* Fourth step: swap columns */
        if (reduction_type == POT_LLL) {
            // Pot-LLL
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
            insert_pos = pot_idx;
        } else {
            if (reduction_type == CLASSIC_LLL) {
                // Standard LLL
                i = (k > low) ? k - 1 : low;
                r_new = R[k][k] * R[k][k] + R[k][i] * R[k][i];
                deep_size = 2;
            } else {
                // Deep insert
                i = low;
                #if BLAS
                    r_new = cblas_ddot(k + 1, R[k], 1, R[k], 1);
                #else
                    for (j = 0, r_new = 0.0; j <= k; ++j) {
                        r_new += R[k][j] * R[k][j];
                    }
                #endif
                deep_size = reduction_type;
            }

            insert_pos = k;
            while (i < k) {
                r_act = delta * R[i][i] * R[i][i];
                //if (delta * R[i][i]*R[i][i] > rhs) {
                if (0.8 * r_act > r_new ||                      // outside of deep-insert slices
                    (
                        (i < deep_size || k - i < deep_size) && // inside of deep-insert slices
                        r_new < r_act
                    )
                   ) {
                    insert_pos = i;
                    break;
                }
                r_new -= R[k][i]*R[k][i];
                i++;
            }
        }

        if (insert_pos < k) {
            if (word_len == WORDLEN_MPZ) {
                swapvl = b[k];
                for (j = k; j > insert_pos; --j) {
                    b[j] = b[j - 1];
                }
                b[insert_pos] = swapvl;
            } else {
                swap = bl[k];
                for (j = k; j > insert_pos; --j) {
                    bl[j] = bl[j - 1];
                }
                bl[insert_pos] = swap;
            }

            // fprintf(stderr, "Swap / rotate %d and %d\n", insert_pos, k);
            k = insert_pos;
        } else {
            k++;
        }
    }
    if (word_len == WORDLEN_MPZ) {
        mpz_clear(musvl);
    }

    if (max_tricols > 1) {
        fprintf(stderr, "Max tricol iterations: %d\n", max_tricols);
    }
    return lowest_pos;

}

void householder_column_inner_hiprec(DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int l, int z, int bit_size) {
    int i;
    DOUBLE w, w_beta, mu;
    DOUBLE eps = 1.0e-15; // 0.0000000001;

    // Apply Householder vectors H[0],..., H[k-1]
    // to R[k] = b[k]:
    //   R[k] -= beta_i * <R[k], H[i]> H[i]
    #if 1
        for (i = 0; i < k; ++i) {
            // w = < R[k], H[i] >
            w = hiprec_dot2(&(R[k][i]), &(H[i][i]), z - i);
            w_beta = w * beta[i];

            #if BLAS
                cblas_daxpy(z - i, -w_beta, &(H[i][i]), 1, &(R[k][i]), 1);
            #else
                for (int j = i; j < z; ++j) {
                    R[k][j] -= w_beta * H[i][j];
                }
            #endif
        }
    #else
        for (i = 0; i < k; ++i) {
            R[i][z - 1] = beta[i] * hiprec_dot2(&(R[k][i]), &(H[i][i]), z - i);
        }
        for (i = 0; i < z; ++i) {
            R[k][i] -= hiprec_dot2_row(&(R[i][z - 1]), z, &(H[i][i]), z, k - i);
        }
    #endif

    // Determine beta and Householder vector H[k]

    // H[k][k] = 1.0;   // Golub, van Loan
    H[k][k] = R[k][k];  // Higham, p.356
    #if BLAS
        cblas_dcopy(z - k - 1, &(R[k][k + 1]), 1, &(H[k][k + 1]), 1);
    #else
        for (i = k + 1; i < z; ++i) {
            H[k][i] = R[k][i];
        }
    #endif

    // Golub, van Loan approach:
    // It ensures that R[k][k] becomes >= 0
    // Higham, p.356: this might cause stability problems
    // DOUBLE sigma, sq;
    // sigma = hiprec_normsq_l2(&(R[k][k + 1]), z - k - 1);
    // if (sigma == 0.0) {
    //     beta[k] = 0.0;
    // } else {
    //     mu = SQRT (R[k][k] * R[k][k] + sigma); // mu = ||R[k]||
    //     if (R[k][k] < -eps) {
    //         H[k][k] = R[k][k] - mu;
    //     } else {
    //         H[k][k] = -sigma / (R[k][k] + mu);
    //     }
    //     sq = H[k][k] * H[k][k];
    //     beta[k] = 2.0 * sq / ( sigma + sq);
    //     for (j = k + 1; j < z; ++j) {
    //         H[k][j] /= H[k][k];
    //     }
    //     H[k][k] = 1.0;
    // }

    // More stable suggestion from Higham:
    mu = hiprec_norm_l2(&(R[k][k]), z - k);
    if (R[k][k] < -eps) {
        mu = -mu;
    }
    H[k][k] += mu;
    beta[k] = 1.0 / (mu * H[k][k]);

    // fprintf(stderr, "H[%d]:", k);
    // for (i = 0; i < z; ++i) {
    //     fprintf(stderr, "%0.20lf ", H[k][i]);
    // }
    // fprintf(stderr, "\n");

    #if 0
    // Check
    if (k >= 0) {
        // Rotate vector R[k] in order to check if
        // H[k] is good enough
        w = hiprec_dot2(&(R[k][k]), &(H[k][k]), z - k);
        w_beta = -w * beta[i];

        for (i = k; i < z; ++i) {
            R[k][i] += w_beta * H[k][i];
        }

        fprintf(stderr, "k: %d, -mu: %0.20lf\n", k, -mu);
        for (i = k; i < z; ++i) {
            fprintf(stderr, "%0.20lf ", R[k][i]);
        }
        fprintf(stderr, "\n");
        // exit(1);
    }
    #endif
    // Apply rotation to R[k][k] only
    // R[k][k] = mu;  // Golub, van Loan
    R[k][k] = -mu;    // Higham
}

int householder_column(mpz_t **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z, int bit_size) {
    int l, j;
    DOUBLE min_val = 0.0;
    DOUBLE min_idx = -1;

    for (l = k; l < s; l++) {
        for (j = 0; j < z; ++j) {
            R[k][j] = (DOUBLE)mpz_get_d(b[l][j]);
        }

        householder_column_inner_hiprec(R, H, beta, k, l, z, bit_size);

        if (l == k || R[k][k] * R[k][k] < min_val) {
            min_val = R[k][k] * R[k][k];
            min_idx = l;
        }
    }
    return min_idx;
}

int householder_column_long(long **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z, int bit_size) {
    int l, j;
    DOUBLE min_val = 0.0;
    DOUBLE min_idx = -1;

    for (l = k; l < s; l++) {
        for (j = 0; j < z; ++j) {
            R[k][j] = (DOUBLE)b[l][j];
        }

        householder_column_inner_hiprec(R, H, beta, k, l, z, bit_size);

        if (l == k || R[k][k] * R[k][k] < min_val) {
            min_val = R[k][k] * R[k][k];
            min_idx = l;
        }
    }
    return min_idx;
}

void size_reduction(mpz_t **b, DOUBLE  **R, mpz_t musvl, double mus, int k, int j, int z) {
    int i;

    for (i = 0; i < z; ++i) {
        mpz_submul(b[k][i], b[j][i], musvl);
    }

    #if BLAS
        cblas_daxpy(j + 1, -mus, R[j], 1, R[k], 1);
    #else
        for (i = 0; i <= j; i++) R[k][i] -= R[j][i] * mus;
    #endif
}

void size_reduction_long(long **b, DOUBLE  **R, long musl, double mus, int k, int j, int z) {
    int i;

    for (i = 0; i < z; ++i) {
        b[k][i] -= musl * b[j][i];
    }

    #if BLAS
        cblas_daxpy(j + 1, -mus, R[j], 1, R[k], 1);
    #else
        for (i = 0; i <= j; i++) R[k][i] -= R[j][i] * mus;
    #endif
}

void check_precision(mpz_t *b, DOUBLE *R, int z, int k) {
    int j;
    mpz_t b_norm;
    DOUBLE r_norm;

    mpz_init(b_norm);
    for (j = 0; j < z; ++j) {
        mpz_addmul(b_norm, b[j], b[j]);
    }
    for (j = 0, r_norm = 0.0; j <= k; ++j) {
        r_norm += R[j] * R[j];
    }
    if (1 || fabs(mpz_get_d(b_norm) - r_norm) > 0.1) {
        fprintf(stderr, "precision check at %d: ", k);
        mpz_out_str(stderr, 10, b_norm);
        fprintf(stderr, " %lf\n", r_norm);
        fflush(stderr);
    }
    mpz_clear(b_norm);
}
