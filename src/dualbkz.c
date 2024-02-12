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
#include "dualbkz.h"

DOUBLE self_dual_bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, DOUBLE p,
                     int (*solutiontest)(lattice_t *lattice, int k)) {
    DOUBLE **R     = lattice->decomp.R;
    DOUBLE *h_beta = lattice->decomp.c;
    // DOUBLE *N      = lattice->decomp.N;
    DOUBLE **H     = lattice->decomp.H;

    DOUBLE r_tt;
    DOUBLE new_cj;
    DOUBLE lD;

    static mpz_t hv;
    int cnt = 0;
    int h, i, last, h_end;
    int start_block, end_block;
    int bit_size = get_bit_size(lattice);

    long *u;

    // Helper arrays for enumerate()
    bkz_enum_t bkz_enum;

    mpz_init(hv);

    last = s - 1; /* |last| points to the last nonzero vector of the lattice.*/
    if (last < 1) {
        printf("BKZ: the number of basis vectors is too small.\n");
        printf("Probably the number of rows is less or equal");
        printf(" to number of columns in the original system\n");
        printf("Maybe you have to increase c0 (the first parameter)!\n");

        mpz_clear(hv);
        return 0.0;
    }

    fprintf(stderr, "\n######### SELFDUAL BKZ ########\n");
    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) {
        u[i] = 0;
    }
    allocate_bkz_enum(&bkz_enum, s);

    while (1) {
        fprintf(stderr, "Start tour #no %d\n", cnt);
        lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, WORDLEN_MPZ, solutiontest);

        fprintf(stderr, "Primal\n");
        //for (start_block = 0; start_block + beta - 1 <= last; ++start_block) {
        //    end_block = start_block + beta - 1;
        for (start_block = 0; start_block < last; ++start_block) {
            end_block = start_block + beta - 1;
            if (end_block > last) end_block = last;

            new_cj = enumerate(lattice, R, u, s, start_block, end_block, delta, p, &bkz_enum);
            h = (start_block - 1 < 0) ? 0 : start_block - 1;
            h_end = (end_block + 1 <= last) ? end_block + 1 : last;

            r_tt = R[start_block][start_block];
            r_tt *= r_tt;
            if (delta * r_tt > new_cj) {
                fprintf(stderr, "primal enumerate successful %d %lf improvement: %lf\n",
                        start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
                fflush(stderr);
                insert_vector(lattice, u, start_block, end_block, z, hv);
                lllH(lattice, R, h_beta, H, h, 0, h_end, z, delta, CLASSIC_LLL, bit_size, WORDLEN_MPZ, solutiontest);
            } else {
                lllH(lattice, R, h_beta, H, h, h, h_end, z, 0.0, CLASSIC_LLL, bit_size, WORDLEN_MPZ, solutiontest);
            }
        }

        break;
        fprintf(stderr, "Dual\n");
        for (start_block = last - beta + 1; start_block > 0; --start_block) {
            end_block = start_block + beta - 1;

            new_cj = dual_enumerate(lattice, R, u, s, start_block, end_block, delta, p, &bkz_enum);
            h = (start_block - 1 < 0) ? 0 : start_block - 1;
            h_end = (end_block + 1 <= last) ? end_block + 1 : last;

            r_tt = 1.0 / R[end_block][end_block];
            r_tt *= r_tt;
            if (delta * r_tt > new_cj) {
                fprintf(stderr, "dual enumerate successful %d %lf improvement: %lf\n",
                        start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
                fflush(stderr);
                dual_insert_vector(lattice, u, start_block, end_block, z, hv);
                lllH(lattice, R, h_beta, H, h, 0, h_end, z, delta, CLASSIC_LLL, bit_size, WORDLEN_MPZ, solutiontest);
            } else {
                lllH(lattice, R, h_beta, H, h, h, h_end, z, 0.0, CLASSIC_LLL, bit_size, WORDLEN_MPZ, solutiontest);
            }
        }

        cnt++;
        if (cnt > 3) break;
    } /* end of |while| */

    lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, WORDLEN_MPZ, solutiontest);

    lD = log_potential(R, s, z);

    #if 1
    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stderr);
    #endif

    free(u);
    free_bkz_enum(&bkz_enum);
    mpz_clear(hv);

    return lD;
}

DOUBLE dual_bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, DOUBLE p,
                int (*solutiontest)(lattice_t *lattice, int k)) {
    DOUBLE **R     = lattice->decomp.R;
    DOUBLE *h_beta = lattice->decomp.c;
    // DOUBLE *N      = lattice->decomp.N;
    DOUBLE **H     = lattice->decomp.H;

    DOUBLE r_tt;
    DOUBLE new_cj, new_cj2;
    DOUBLE lD;

    static mpz_t hv;
    int cnt;
    int h, i, last, h_end;
    int start_block, end_block;
    int bit_size = get_bit_size(lattice);

    long *u;
    // Helper arrays for enumerate()
    bkz_enum_t bkz_enum;

    mpz_init(hv);

    last = s - 1; /* |last| points to the last nonzero vector of the lattice.*/
    if (last < 1) {
        printf("BKZ: the number of basis vectors is too small.\n");
        printf("Probably the number of rows is less or equal");
        printf(" to number of columns in the original system\n");
        printf("Maybe you have to increase c0 (the first parameter)!\n");

        mpz_clear(hv);
        return 0.0;
    }

    fprintf(stderr, "\n######### DUAL BKZ ########\n");
    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) {
        u[i] = 0;
    }
    allocate_bkz_enum(&bkz_enum, s);

    //decomp_alloc(&R, &h_beta, &N, &H, s, z);
    lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, WORDLEN_MPZ, solutiontest);

    cnt = -1;
    end_block = last + 1;
    //start_block = 0;
    while (cnt < last) {
        end_block--;
        if (end_block < 2) {
            break;
        }
        h_end = (end_block < last) ? end_block + 1 : last;

        start_block = end_block - beta + 1;
        start_block = (start_block >= 0) ? start_block : 0;

        new_cj = dual_enumerate(lattice, R, u, s, start_block, end_block, delta, p, &bkz_enum);
        h = (start_block - 1 < 0) ? 0 : start_block - 1;

        r_tt = 1.0 / R[end_block][end_block];
        r_tt *= r_tt;
        if (delta * r_tt > new_cj) {
            fprintf(stderr, "dual enumerate successful %d %lf improvement: %lf\n",
                    start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
            fflush(stderr);

            /* successful enumeration */
            dual_insert_vector(lattice, u, start_block, end_block, z, hv);
            lllH(lattice, R, h_beta, H, h, h, h_end, z, 0.0, CLASSIC_LLL, bit_size, WORDLEN_MPZ, solutiontest);
            i = end_block;

            new_cj2 = 1.0 / (R[i][i] * R[i][i]);
            if (FALSE && fabs(new_cj2 - new_cj) > EPSILON) {
                fprintf(stderr, "???????????????? We have a problem at %d: %lf %lf\n", i, new_cj2, new_cj);
                fflush(stderr);
                exit(EXIT_ERR_INPUT);
            }
            lllH(lattice, R, h_beta, H, h, 0, h_end, z, delta, CLASSIC_LLL, bit_size, WORDLEN_MPZ, solutiontest);
            //cnt = -1;
            cnt++;
        } else {
            cnt++;
        }
    } /* end of |while| */

    lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, WORDLEN_MPZ, solutiontest);

    lD = log_potential(R, s, z);

    #if FALSE
    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stderr);
    #endif

    free(u);
    free_bkz_enum(&bkz_enum);

    mpz_clear(hv);

    return lD;
}

DOUBLE dual_enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s,
                      int start_block, int end_block, DOUBLE improve_by, DOUBLE p,
                      bkz_enum_t *bkz_enum) {

    DOUBLE *y, *c, *a;
    long *delta, *d, *v;
    DOUBLE *u_loc;
    DOUBLE c_min;

    int i, j;
    int t, t_min;
    int found_improvement = 0;

    int len, k;
    double alpha, radius;
    int SCHNITT = 2000;

    c = bkz_enum->c;
    y = bkz_enum->y;
    d = bkz_enum->d;
    v = bkz_enum->v;
    delta = bkz_enum->delta;
    u_loc = bkz_enum->u_loc;
    a = (DOUBLE*)calloc(s + 1, sizeof(DOUBLE));

    len = end_block + 1 - start_block;
    for (i = start_block; i <= end_block; i++) {
        c[i] = y[i] = a[i] = 0.0;
        u_loc[i] = 0.0;
        v[i] = delta[i] = 0;
        d[i] = 1;
    }
    c_min = 1.0 / R[end_block][end_block];
    c_min *= c_min;
    c_min *= improve_by;

    for (t_min = end_block - 1; t_min >= start_block; t_min--) {
        t = t_min;
        u_loc[t] = 1.0;
        len = end_block + 1 - t_min;

        //fprintf(stderr, "LOOP %d %d %d\n", t_min, end_block, s+1);
        while (t >= t_min) {
            handle_signals(lattice, R);

            a[t] = u_loc[t] - y[t];
            //    c[t] = c[t - 1] + (a[t]/ R[t][t])^2
            c[t] = a[t] / R[t][t];
            c[t] *= c[t];
            if (t > t_min) {
                c[t] += c[t - 1];
            }

            if (len <= SCHNITT) {
                alpha = 1.0;
            } else {
                k = t - start_block + 1;
                if (k > 1 * len / 3) {
                    alpha = 1.0;
                } else {
                    //alpha = p;
                    alpha = 3 * p * k / len;
                }
                alpha = (alpha < 1.0) ? alpha : 1.0;
            }
            radius = alpha * c_min;

            if (c[t] < radius - EPSILON) {
                if (t < end_block) {
                    // forward
                    t++;

    #if BLAS
                    y[t] = cblas_ddot(t - t_min, &(a[t_min]), 1, &(R[t][t_min]), 1);
    #else
                    for (j = t_min, y[t] = 0.0; j < t; j++) {
                        y[t] += a[j] * R[t][j];
                    }
    #endif
                    y[t] /= R[t][t];

                    u_loc[t] = v[t] = (long)(ROUND(y[t]));
                    delta[t] = 0;
                    d[t] = (v[t] > y[t]) ? -1 : 1;

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
                t--;
                if (t < t_min) {
                    break;
                }
            }
            // next
            if (t == t_min) {
                u_loc[t] += 1;
                //fprintf(stderr, "%d %ld %ld u=%lf v=%ld, %ld\n", t, d[t], delta[t], u_loc[t], v[t], v[t] + delta[t]);
            } else {
                if (t > t_min) delta[t] *= -1.0;
                if (delta[t] * d[t] >= 0) delta[t] += d[t];
                u_loc[t] = v[t] + delta[t];
            }
        }
    }

    free(a);

    if (!found_improvement) {
        c_min = 1.0 / R[end_block][end_block];
        c_min *= c_min;
    }
    return (c_min);
}

void dual_insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv) {
    mpz_t **b = lattice->basis;
    mpz_t *swapvl;
    int i, j;
    // int g;
    // long q ui;

    /* build new basis */
    for (j = 0; j < z; j++) {
        mpz_set_si(lattice->swap[j], 0);
    }

    // Store new linear combination in lattice->swap
    for (i = start; i <= end; i++) {
        if (u[i] != 0) {
            for (j = 0; j < z; j++) {
                // There is no mpz_addmul_si
                if (u[i] > 0) {
                    mpz_addmul_ui(lattice->swap[j], b[i][j], u[i]);
                } else {
                    mpz_submul_ui(lattice->swap[j], b[i][j], -u[i]);
                }
            }
        }
    }

    swapvl = b[lattice->num_cols];
    for (i = lattice->num_cols; i > start; i--) {
        b[i] = b[i - 1];
    }
    b[start] = lattice->swap;
    lattice->swap = swapvl;
    lattice->num_cols++;
}
