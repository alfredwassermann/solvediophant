#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "const.h"
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

/**
 *  Lattice basis reduction algorithms
 */
int lllH(lattice_t *lattice, DOUBLE **R, DOUBLE *beta, DOUBLE **H,
            int start, int low, int up, int z,
            DOUBLE delta, int reduction_type,
            int bit_size,
            int (*solutiontest)(lattice_t *lattice, int k)) {

    coeff_t **b = lattice->basis;
    int i, j, k;
    int cnt_tricol;
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
    int deep_size;
    coeff_t *swapvl;

    #if VERBOSE > 1
        int counter = 0;
    #endif

    lattice->work_on_long = FALSE;

    eta = ETACONST;
    if (bit_size > 100) {
        theta = 0.50;
        eta = 0.52;
    } else if (bit_size > 55) {
        eta = 0.52;
        theta = 0.3;
    } else if (bit_size > 30) {
        theta = 0.04;
    } else {
        theta = 0.0;
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

    mpz_init(musvl);
    mpz_init(hv);
    //mpz_init(sum_mu);

    /* Test for trivial cases. */
    if ((z <= 1) || (up <= 1)) {
        fprintf(stderr, "Wrong dimensions in lllH\n");
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
        #if VERBOSE > 1
            if (counter > 0 && (counter % 10000) == 0) {
                fprintf(stderr, "LLL: %d k:%d\n", counter, k);
                fflush(stderr);
            }
            counter++;
        #endif
        handle_signals(lattice, R);

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


        cnt_tricol = 0;
    start_tricol:
        /* Recompute column k of R */
        i = householder_column(b, R, H, beta, k, k + 1, z, bit_size);

        /* size reduction of $b_k$ */
        mu_all_zero = TRUE;
        for (j = k - 1; j >= low; j--) {
            /* Subtract suitable multiple of $b_j$ from $b_k$. */
            if (fabs(R[k][j]) > eta * fabs(R[j][j]) + theta * fabs(R[k][k])) {
                mus = ROUND(R[k][j] / R[j][j]);
                mpz_set_d(musvl, mus);
                mu_all_zero = FALSE;

                if (cnt_tricol > 1000) {
                    fprintf(stderr, "Possible tricol error: %d: eta=%0.2lf, theta=%0.2lf, %0.2lf, %lf %lf %lf\n\t %lf\n", 
                        j, eta, theta, mus,
                        R[k][j], R[j][j], R[k][k],
                        eta * fabs(R[j][j]) + theta * fabs(R[k][k]));
                    //exit(1);
                }
                /* set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$ */
                size_reduction(b, R, musvl, mus, k, j);
                (*solutiontest)(lattice, k);
            } else {}
        }

        if (cnt_tricol > 0 && cnt_tricol % 1000 == 0) {
            fprintf(stderr, "tricol %d at %d\n", cnt_tricol, k);
            fflush(stderr);
        }
        cnt_tricol++;

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
            matrix and restart lllH with s = s-1.
        */
        #if BLAS
            norm = cblas_dnrm2(k + 1, R[k], 1);
        #else
            for (j = 0, norm = 0.0; j <= k; ++j) {
                norm += R[k][j] * R[k][j];
            }
            norm = SQRT(norm);
        #endif

        if (norm != norm || norm < 0.5) {  // nan or < 0.5
            //print_lattice(lattice);
            swapvl = b[k];
            //for (i = k + 1; i < up; i++) {
            for (i = k + 1; i < lattice->num_cols; i++) {
                b[i-1] = b[i];
            }
            b[lattice->num_cols] = swapvl;

            up--;
            lattice->num_cols--;

            k = 0;
            //print_lattice(lattice);
            fprintf(stderr, "Zero vector at %d\n", k);
            fflush(stderr);
            continue;
        }

        // If delta == 0, only size reduction is done
        if (delta == 0.0) {
            k++;
            continue;
        }

        /* fourth step: swap columns */
        if (reduction_type != POT_LLL) {
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
                 if (0.8 * r_act > r_new ||
                     ((i < deep_size || k - i < deep_size) &&
                      r_act > r_new)) {
                    insert_pos = i;
                    break;
                 }
                 r_new -= R[k][i]*R[k][i];
                 i++;
            }
        } else {
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
        }

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
                norm += R[k][j] * R[k][j];
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
        for (i = 0; i < j; i++) mu[k][i] += mu[j][i];
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
        for (i = 0; i < j; i++) mu[k][i] -= mu[j][i] * mus;
    #endif

    }
}

int lllH_long(lattice_t *lattice, DOUBLE **R, DOUBLE *beta, DOUBLE **H,
            int start, int low, int up, int z,
            DOUBLE delta, int reduction_type,
            int bit_size,
            int (*solutiontest)(lattice_t *lattice, int k)) {

    long **b = lattice->basis_long;
    int i, j, k;
    int cnt_tricol;
    DOUBLE norm;
    int mu_all_zero;

    DOUBLE theta, eta;

    DOUBLE mus;
    long musl;

    DOUBLE r_new, r_act;
    DOUBLE pot, pot_max;
    int pot_idx;
    int insert_pos, lowest_pos;
    int deep_size;
    long *swap;

    #if VERBOSE > 1
        int counter = 0;
    #endif

    lattice->work_on_long = TRUE;

    eta = ETACONST;
    if (bit_size > 100) {
        theta = 0.50;
        eta = 0.52;
    } else if (bit_size > 55) {
        theta = 0.1;
    } else if (bit_size > 30) {
        theta = 0.04;
    } else {
        theta = 0.0;
    }

    #if VERBOSE > 1
        fprintf(stderr, "LLLH_long: ");
        if (reduction_type == POT_LLL) {
            fprintf(stderr, "use PotLLL");
        } else if (reduction_type == CLASSIC_LLL) {
            fprintf(stderr, "use classic LLL");
        } else {
            fprintf(stderr, "use deepinsert %d", reduction_type);
        }
        fprintf(stderr, ", delta=%0.3lf", delta);
        fprintf(stderr, ", eta=%0.3lf", eta);
        fprintf(stderr, ", max bits: %d\n", Pobit_size);
    #endif

    /* Test for trivial cases. */
    if ((z <= 1) || (up <= 1)) {
        fprintf(stderr, "Wrong dimensions in LLLH_long\n");
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
        #if VERBOSE > 1
            if (counter > 0 && (counter % 10000) == 0) {
                fprintf(stderr, "LLL: %d k:%d\n", counter, k);
                fflush(stderr);
            }
            counter++;
        #endif
        handle_signals(lattice, R);

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


        cnt_tricol = 0;
    start_tricol:
        /* Recompute column k of R */
        i = householder_column_long(b, R, H, beta, k, k + 1, z, bit_size);

        /* size reduction of $b_k$ */
        mu_all_zero = TRUE;
        for (j = k - 1; j >= low; j--) {
            /* Subtract suitable multiple of $b_j$ from $b_k$. */
            if (fabs(R[k][j]) > eta * fabs(R[j][j]) + theta * fabs(R[k][k])) {
                mus = ROUND(R[k][j] / R[j][j]);
                musl = (long)mus;
                mu_all_zero = FALSE;

                if (cnt_tricol > 1000) {
                    fprintf(stderr, "%d: eta=%0.2lf, theta=%0.2lf, %0.2lf, %lf %lf %lf\n\t %lf\n", j, eta, theta, mus,
                        fabs(R[k][j]), fabs(R[j][j]), fabs(R[k][k]),
                        eta * fabs(R[j][j]) + theta * fabs(R[k][k]));
                }
                /* set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$ */
                size_reduction_long(b, R, musl, mus, k, j, z);
                (*solutiontest)(lattice, k);
            } else {}
        }

        if (cnt_tricol > 0 && cnt_tricol % 1000 == 0) {
            fprintf(stderr, "tricol %d at %d\n", cnt_tricol, k);
            fflush(stderr);
        }
        cnt_tricol++;

        if (!mu_all_zero) {
            goto start_tricol;
        }

        /*
            Before going to step 4 we test if $b_k$ is linear dependent.
            If we find a linear dependent vector $b_k$,
            we shift b_k to the last column of the
            matrix and restart lllH with s = s-1.
        */
        #if BLAS
            norm = cblas_dnrm2(k + 1, R[k], 1);
        #else
            for (j = 0, norm = 0.0; j <= k; ++j) {
                norm += R[k][j] * R[k][j];
            }
            norm = SQRT(norm);
        #endif

        if (norm != norm || norm < 0.5) {  // nan or < 0.5
            //print_lattice(lattice);
            swap = b[k];
            //for (i = k + 1; i < up; i++) {
            for (i = k + 1; i < lattice->num_cols; i++) {
                b[i-1] = b[i];
            }
            b[lattice->num_cols] = swap;

            up--;
            lattice->num_cols--;

            k = 0;
            //print_lattice(lattice);
            fprintf(stderr, "Zero vector at %d\n", k);
            fflush(stderr);
            continue;
        }

        // If delta == 0, only size reduction is done
        if (delta == 0.0) {
            k++;
            continue;
        }

        /* fourth step: swap columns */
        if (reduction_type != POT_LLL) {
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
                 if (0.8 * r_act > r_new ||
                     ((i < deep_size || k - i < deep_size) &&
                      r_act > r_new)) {
                    insert_pos = i;
                    break;
                 }
                 r_new -= R[k][i]*R[k][i];
                 i++;
            }
        } else {
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
        }

        if (insert_pos < k) {
            swap = b[k];
            for (j = k; j > insert_pos; --j) {
                b[j] = b[j - 1];
            }
            b[insert_pos] = swap;

            //fprintf(stderr, "INSERT %d at %d\n", k, insert_pos);
            k = insert_pos;
        } else {
            k++;
        }
    }

    return lowest_pos;

}

int householder_column_long(long **b, DOUBLE **R, DOUBLE **H, DOUBLE *beta, int k, int s, int z, int bit_size) {
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
            R[k][j] = b[l][j];
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
                norm += R[k][j] * R[k][j];
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

void size_reduction_long(long **b, DOUBLE  **mu, long musl, double mus, int k, int j, int z) {
    int i;

    for (i = 0; i < z; ++i) {
        b[k][i] -= musl * b[j][i];
    }
    #if BLAS
        cblas_daxpy(j, -mus, mu[j], 1, mu[k], 1);
    #else
        for (i = 0; i < j; i++) mu[k][i] -= mu[j][i] * mus;
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
