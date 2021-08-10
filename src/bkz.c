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
 * Blockwise Korkine Zolotareff reduction
 */
DOUBLE bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, DOUBLE p,
            int enum_type, int max_tours,
            int (*solutiontest)(lattice_t *lattice, int k),
            int (*solutiontest_long)(lattice_t *lattice, int k)) {

    DOUBLE **R     = lattice->decomp.R;
    DOUBLE *h_beta = lattice->decomp.c;
    DOUBLE **H     = lattice->decomp.H;
    DOUBLE r_tt;
    DOUBLE new_cj;
    DOUBLE lD;

    static mpz_t hv;
    int cnt;
    int tour_cnt, enum_cnt;
    int h, i, last;
    int start_block, end_block;
    int bit_size = get_bit_size(lattice);

    long *u;
    int bit_size_threshold = 32;

    // Helper arrays for enumerate()
    bkz_enum_t bkz_enum;

    int j;
    for (i = 0; i < lattice->num_cols; i++) {
        for (j = 0; j < lattice->num_rows; j++) {
            H[i][j] = 0.0;
        }
        for (j = 0; j < lattice->num_cols; j++) {
            R[i][j] = 0.0;
        }
        h_beta[i] = 0.0;
    }
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

    if (enum_type == ENUM_LDS_FULL) {
        fprintf(stderr, "\n######### BKZ-LDS %d ########\n", beta);
    } else if (enum_type == ENUM_LDS_FULL2) {
        fprintf(stderr, "\n######### BKZ-LDS2 %d ########\n", beta);
    } else {
        fprintf(stderr, "\n######### BKZ %d ########\n", beta);
    }

    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) {
        u[i] = 0;
    }

    if (bit_size < bit_size_threshold) {
        copy_lattice_to_long(lattice);
        lllH_long(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest_long);
    } else {
        lllH     (lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest);
    }

    allocate_bkz_enum(&bkz_enum, s);

    start_block = cnt = -1;
    tour_cnt = 0;
    enum_cnt = 0;

    // print_gsa(R, s);

    while (cnt < last && tour_cnt < max_tours) {
        start_block++;

        if (start_block == last) {
            start_block = 0;
            tour_cnt++;
        }
        if (start_block == 0) {
            //lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size);
        }

        end_block = start_block + beta - 1;
        end_block = (end_block < last) ? end_block : last;

        enum_cnt++;
        if (enum_type != ENUM_BLOCK) {
            //fprintf(stderr, "LDS %d\n", tour_cnt);
            // fprintf(stderr, "lds at %d, \tstart at %d (%d)\n", cnt, start_block, tour_cnt);
            if (enum_type == ENUM_LDS_FULL || enum_type == ENUM_LDS_FULL2) {
                end_block = last;
            }
            new_cj = lds_enumerate(lattice, R, u, s, start_block, end_block, delta, p, &bkz_enum);
        } else {
            //fprintf(stderr, "BLOCK %d\n", tour_cnt);
            // fprintf(stderr, "block at %d, \tstart at %d (%d)\n", cnt, start_block, tour_cnt);
            new_cj = enumerate(lattice, R, u, s, start_block, end_block, delta, p, &bkz_enum);
        }

        h = (end_block + 1 < last) ? end_block + 1 : last;

        r_tt = R[start_block][start_block];
        r_tt *= r_tt;
        if (delta * r_tt > new_cj) {
            #if FALSE
            if (beta > 20) {
                fprintf(stderr, "enum %d successful %d %lf improvement: %lf\n",
                                    beta, start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
                fflush(stderr);
            }
            #endif

            /* successful enumeration */
            if (bit_size < bit_size_threshold) {
                // fprintf(stderr, "long insert\n");
                insert_vector_long(lattice, u, start_block, end_block, z);
            } else {
                // fprintf(stderr, "mpz insert\n");
                insert_vector     (lattice, u, start_block, end_block, z, hv);
            }

            /*
            i = householder_column_long(lattice->basis_long, R, H, h_beta, start_block, start_block + 1, z, bit_size);
            new_cj2 = R[i][i] * R[i][i];
            if (fabs(new_cj2 - new_cj) > EPSILON) {
                fprintf(stderr, "???????????????? We have a problem: %lf %lf\n", new_cj2, new_cj);
                fflush(stdout);
            }
            */

            if (bit_size < bit_size_threshold) {
                lllH_long(lattice, R, h_beta, H, start_block - 1, 0, /*lattice->num_cols */ h + 1, z, delta, POT_LLL, bit_size, solutiontest_long);
            } else {
                lllH     (lattice, R, h_beta, H, start_block - 1, 0, /*lattice->num_cols */ h + 1, z, delta, POT_LLL, bit_size, solutiontest);
            }
            //lllH(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, 0.0, CLASSIC_LLL, bit_size, solutiontest);
            //lattice->num_cols--;

            //start_block = lllH(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, delta, CLASSIC_LLL, bit_size);
            //fprintf(stderr, "%d\n", start_block);

            if (enum_type == ENUM_LDS_FULL2) {
                start_block--;
                if (enum_cnt > 30000) {
                    break;
                }
            }

            // print_gsa(R, s);
            cnt = -1;

        } else {
            // fprintf(stderr, "enumerate: no improvement %d\n", cnt);
            // fflush(stderr);

            if (h > 0) {
                if (bit_size < bit_size_threshold) {
                    lllH_long(lattice, R, h_beta, H, h - 1, h - 1, h + 1, z, 0.0, CLASSIC_LLL, bit_size, solutiontest_long);
                } else {
                    lllH     (lattice, R, h_beta, H, h - 1, h - 1, h + 1, z, 0.0, CLASSIC_LLL, bit_size, solutiontest);
                }
            }
            cnt++;
        }

    } /* end of |while| */

    if (bit_size < bit_size_threshold) {
        copy_lattice_to_mpz(lattice);
    }

    lD = log_potential(R, s, z);
    free_bkz_enum(&bkz_enum);

    #if 0
    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stderr);
    #endif

    free(u);
    mpz_clear(hv);

    return lD;
}

/**
 * @brief Insert a linear combination of lattice vectors at position `start` and
 *        shift all previous lattice vectors of the whole lattice starting form this position by one.
 *        Increase lattice->num_cols.
 *        Use GMP to compute the linear combination.
 *
 * @param {lattice_t} lattice
 * @param {long*} u Vector containing the coefficients of the linear combination
 * @param {int} start Position in the array of lattice vectors where the linear combination starts, and where the
 *              vector will be inserted.
 * @param {int} end Last position (included) of the linear combination.
 * @param {int} z Length of the basis vectors
 * @param {mpz_t} hv Multiprecision helper variable
 */
void insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv) {
    coeff_t **b = lattice->basis;
    coeff_t *swapvl;
    int i, j, g;
    long q, ui;

    /* build new basis */
    for (j = 1; j <= z; j++) {
        mpz_set_si(lattice->swap[j].c, 0);
    }

    // Store new linear combination in lattice->swap
    for (i = start; i <= end; i++) {
        if (u[i] != 0) for (j = 1; j <= z; j++) {
            if (u[i] > 0) {
                mpz_addmul_ui(lattice->swap[j].c, b[i][j].c, u[i]);
            } else {
                mpz_submul_ui(lattice->swap[j].c, b[i][j].c, -u[i]);
            }
        }
    }
    coeffinit(lattice->swap, z);

    #if 1
        swapvl = b[lattice->num_cols];
        for (i = lattice->num_cols; i > start; i--) {
            b[i] = b[i - 1];
        }
        b[start] = lattice->swap;
        lattice->swap = swapvl;
        lattice->num_cols++;
    #else
        g = end;
        while (u[g] == 0) g--;
        i = g - 1;
        while (labs(u[g]) > 1) {
            while (u[i] == 0) i--;
            q = (long)ROUND((1.0 * u[g]) / u[i]);
            ui = u[i];
            u[i] = u[g] - q*u[i];
            u[g] = ui;

            // (b[g], b[i]) = (b[g] * q + b[i], b[g])
            for (j = 1; j <= z; j++) {
                mpz_set(hv, b[g][j].c);
                mpz_mul_si(b[g][j].c, b[g][j].c, (long)q);
                mpz_add(b[g][j].c, b[g][j].c, b[i][j].c);
                mpz_set(b[i][j].c, hv);
            }
            coeffinit(b[g], z);
            coeffinit(b[i], z);
        }

        // (b[start], b[start+1], ... , b[g]) -> (b[g], b[start], ... , b[g-1])
        swapvl = b[g];
        for (i = g; i > start; i--) {
            b[i] = b[i - 1];
        }
        b[start] = lattice->swap;
        coeffinit(b[start], z);

        lattice->swap = swapvl;
        for (j = 1; j <= z; j++)
            mpz_set_si(lattice->swap[j].c, 0);
        coeffinit(lattice->swap, z);

        #if 0
        for (j = 0; j < z; j++) {
            mpz_out_str(stderr, 10, get_entry(lattice->basis, start, j));
            fprintf(stderr," ");
        }
        fprintf(stderr, "\n");
        fflush(stderr);
        #endif
    #endif
}

/**
 * @brief Insert a linear combination of lattice vectors at position `start` and
 *        shift all previous lattice vectors of the whole lattice starting form this position by one.
 *        Increase lattice->num_cols.
 *
 * @param {lattice_t} lattice
 * @param {long*} u Vector containing the coefficients of the linear combination
 * @param {int} start Position in the array of lattice vectors where the linear combination starts, and where the
 *              vector will be inserted.
 * @param {int} end Last position (included) of the linear combination.
 * @param {int} z Length of the basis vectors
 */
void insert_vector_long(lattice_t *lattice, long *u, int start, int end, int z) {
    long **b = lattice->basis_long;
    long *swap;
    int i, j, g;
    long q, ui;
    long hv;

    // Store new linear combination in lattice->swap_long
    for (j = 0; j < z; j++) {
        lattice->swap_long[j] = 0;
    }

    for (i = start; i <= end; i++) {
        if (u[i] != 0) for (j = 0; j < z; j++) {
            lattice->swap_long[j] += b[i][j] * u[i];
        }
    }

    #if 1
        swap = b[lattice->num_cols];
        for (i = lattice->num_cols; i > start; i--) {
            b[i] = b[i - 1];
        }
        b[start] = lattice->swap_long;
        lattice->swap_long = swap;

        lattice->num_cols++;
    #else
        g = end;
        while (u[g] == 0) g--;
        i = g - 1;
        while (labs(u[g]) > 1) {
            while (u[i] == 0) i--;
            q = (long)ROUND((1.0 * u[g]) / u[i]);
            ui = u[i];
            u[i] = u[g] - q*u[i];
            u[g] = ui;

            // (b[g], b[i]) = (b[g] * q + b[i], b[g])
            for (j = 0; j < z; j++) {
                hv = b[g][j];
                b[g][j] = b[g][j] * q + b[i][j];
                b[i][j] = hv;
            }
        }

        // (b[start], b[start+1], ... , b[g]) -> (b[g], b[start], ... , b[g-1])
        swap = b[g];
        for (i = g; i > start; i--) {
            b[i] = b[i - 1];
        }
        b[start] = lattice->swap_long;

        lattice->swap_long = swap;
        for (j = 0; j < z; j++) {
            lattice->swap_long[j] = 0;
        }

    #endif
}

void allocate_bkz_enum(bkz_enum_t *bkz_enum, int s) {
    bkz_enum->c = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    bkz_enum->y = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    bkz_enum->d = (long*)calloc(s+1,sizeof(long));
    bkz_enum->v = (long*)calloc(s+1,sizeof(long));
    bkz_enum->delta = (long*)calloc(s+1,sizeof(long));
    bkz_enum->u_loc = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
}

void free_bkz_enum(bkz_enum_t *bkz_enum) {
    free(bkz_enum->c);
    free(bkz_enum->y);
    free(bkz_enum->d);
    free(bkz_enum->v);
    free(bkz_enum->delta);
    free(bkz_enum->u_loc);
}

/**
 * Pruned Gauss-Enumeration.
 */
DOUBLE enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s,
                    int start_block, int end_block, DOUBLE improve_by, DOUBLE p,
                    bkz_enum_t *bkz_enum) {

    DOUBLE *y, *c;
    long *delta, *d, *v;
    DOUBLE *u_loc;
    DOUBLE c_min;

    int i, j;
    int t, t_max;
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

    for (i = start_block; i <= end_block + 1; i++) {
        c[i] = y[i] = 0.0;
        u_loc[i] = 0.0;
        v[i] = delta[i] = 0;
        d[i] = 1;
    }

    if (end_block - start_block <= SCHNITT) {
        c_min = set_prune_const(R, start_block, end_block + 1, PRUNE_NO, 1.0);
    } else {
        //Hoerners version of the Gaussian volume heuristics.
        //hoerner(R, start_block, end_block + 1, p, eta);

        // New heuristic
        c_min = set_prune_const(R, start_block, end_block + 1, PRUNE_BKZ, p);
    }
    c_min *= improve_by;

    //for (t_max = start_block + 1; t_max <= end_block; t_max ++) {
    t = t_max = start_block + 1;
    for (i = start_block; i <= end_block + 1; i++) {
        c[i] = y[i] = 0.0;
        u_loc[i] = 0.0;
        v[i] = delta[i] = 0;
        d[i] = 1;
    }
    // t = t_max;
    u_loc[t] = 1.0;
    len = t_max + 1 - start_block;
    while (t <= end_block) {
        handle_signals(lattice, R);
        /*
            c[t] = ((u_loc[t] + y[t]) * R[t][t])^2 + c[t + 1]
        */
        c[t] = (u_loc[t] + y[t]) * R[t][t];
        c[t] *= c[t];
        c[t] += c[t + 1];
        if (len <= SCHNITT) {
            alpha = 1.0;
        } else {
            k = (end_block + 1 - t);
            if (k > 6 * len / 9) {
                alpha = 1.0;
            } else {
                //alpha = p;
                alpha = 3 * p * k / len;
            }
            alpha = (alpha < 1.0) ? alpha : 1.0;
        }

        //fprintf(stderr, "%d %d %d %lf\n", start_block, t, end_block, alpha);

        radius = alpha * c_min;

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
                    // fprintf(stderr, "%ld ", u[i]);
                }
                // fprintf(stderr, "\n");
                found_improvement = 1;

                // Immediate return after the first improvement,
                // return c_min;
            }
        } else {
            // back
            t++;
            if (t > t_max) {
                t_max = t;
                //break;
            }
        }
        // next
        if (t == t_max) {
            //fprintf(stderr, "%d %ld %ld u=%lf v=%ld, %ld\n",
            //      t, d[t], delta[t], u_loc[t], v[t], v[t] + delta[t]);
            u_loc[t] += 1;
        } else {
            if (t < t_max) delta[t] *= -1.0;
            if (delta[t] * d[t] >= 0) delta[t] += d[t];
            u_loc[t] = v[t] + delta[t];
        }
    }
    //}

    //free(lambda_min);

    if (!found_improvement) {
        c_min = R[start_block][start_block];
        c_min *= c_min;
    }
    return (c_min);
}

DOUBLE lds_enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s,
                    int start_block, int end_block, DOUBLE improve_by, DOUBLE p,
                    bkz_enum_t *bkz_enum) {

    DOUBLE *y, *c;
    long *delta, *d, *v;
    DOUBLE *u_loc;

    DOUBLE c_min;

    int i, j;
    int t, t_max;
    int found_improvement = 0;

    int len, k;
    double alpha, radius;
    int SCHNITT = 2000;

    // --------
    int *lds_k;
    int lds_k_start;
    int lds_k_max;

    lds_k_max = 2;
    lds_k = (int*)malloc((s + 1) * sizeof(int));

    // --------

    if (end_block - start_block < 30) {
        lds_k_max = 9;
    } else if (end_block - start_block < 60) {
        lds_k_max = 7;
    } else if (end_block - start_block < 100) {
        lds_k_max = 4;
    }

    c = bkz_enum->c;
    y = bkz_enum->y;
    d = bkz_enum->d;
    v = bkz_enum->v;
    delta = bkz_enum->delta;
    u_loc = bkz_enum->u_loc;

    for (i = start_block; i <= end_block + 1; i++) {
        c[i] = y[i] = 0.0;
        u_loc[i] = 0.0;
        v[i] = delta[i] = 0;
        d[i] = 1;
    }

    if (end_block - start_block <= SCHNITT) {
        c_min = set_prune_const(R, start_block, end_block + 1, PRUNE_NO, 1.0);
    } else {
        //Hoerners version of the Gaussian volume heuristics.
        //hoerner(R, start_block, end_block + 1, p, eta);

        // New heuristic
        c_min = set_prune_const(R, start_block, end_block + 1, PRUNE_BKZ, p);
    }
    c_min *= improve_by;

    //for (t_max = start_block + 1; t_max <= end_block; t_max ++) {
    for (t_max = end_block; t_max > start_block; t_max--) {

        //t = t_max = start_block + 1;
        for (i = start_block; i <= end_block + 1; i++) {
            c[i] = y[i] = 0.0;
            u_loc[i] = 0.0;
            v[i] = delta[i] = 0;
            d[i] = 1;
        }

        for (lds_k_start = 1; lds_k_start < lds_k_max; lds_k_start++) {
            // fprintf(stderr, "LDS %d\n", lds_k_start);

            t = t_max;
            u_loc[t] = 1.0;
            c[t] = y[t] = 0.0;
            v[t] = delta[t] = 0;
            d[t] = 1;

            if (t < 0 || t > s) {
                fprintf(stderr, "Y%d\n", t);
            }

            lds_k[t] = lds_k_start;

            len = t_max + 1 - start_block;
            while (t <= t_max) {
                //fprintf(stderr, "t=%d, start=%d, end=%d t_max=%d lds:%d\n", t, start_block, end_block, t_max, lds_k[t]);
                handle_signals(lattice, R);
                /*
                    c[t] = ((u_loc[t] + y[t]) * R[t][t])^2 + c[t + 1]
                */
                c[t] = (u_loc[t] + y[t]) * R[t][t];
                c[t] *= c[t];
                c[t] += c[t + 1];

                alpha = 1.0;

                //fprintf(stderr, "%d %d %d %lf\n", start_block, t, end_block, alpha);
                radius = alpha * c_min;
                if (c[t] < radius - EPSILON) {
                    if (t > start_block) {
                        // forward
                        t--;

            if (t < 0 || t+1 > s) {
                fprintf(stderr, "Z%d\n", t);
            }

                        lds_k[t] = lds_k[t + 1];
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
                        }
                        #if FALSE
                        fprintf(stderr, "i ");
                        for (i = start_block; i <= end_block; i++) {
                            fprintf(stderr, "%ld ", u[i]);
                        }
                        fprintf(stderr, "\n");
                        #endif
                        found_improvement = 1;

                        // Immediate return after the first improvement,
                        // return c_min;
                    }

                } else {
                    // back
lds_back:
                    // fprintf(stderr, "back %d\n", t);
                    t++;
                    if (t > t_max) {
                        break;
                    }
                }
                // next
                // fprintf(stderr, "t:%d lds:%d\n", t, lds_k[t]);
                if (t == t_max) {
                    //fprintf(stderr, "%d %ld %ld u=%lf v=%ld, %ld\n",
                    //      t, d[t], delta[t], u_loc[t], v[t], v[t] + delta[t]);
                    u_loc[t] += 1;
                } else {
                    if (t < t_max) delta[t] *= -1.0;
                    if (delta[t] * d[t] >= 0) delta[t] += d[t];
                    u_loc[t] = v[t] + delta[t];
                }
                    if (t < 0 || t > s) {
                        fprintf(stderr, "W%d\n", t);
                    }
                if (lds_k[t] == 0) {
                    // fprintf(stderr, "GOTO\n");
                    goto lds_back;
                } else {
                    if (t < 0 || t > s) {
                        fprintf(stderr, "Z%d\n", t);
                    }
                    lds_k[t]--;
                }
            }
        }
    }

    if (!found_improvement) {
        c_min = R[start_block][start_block];
        c_min *= c_min;
    }

    free(lds_k);

    return (c_min);
}

/*
void set_linear_pruning_const(DOUBLE *alpha, len, ) {
    int i;
    if (len <= SCHNITT) {
        alpha = 1.0;
    } else {
        k = (end_block + 1 - t);
        if (k > 6 * len / 9) {
            alpha = 1.0;
        } else {
            //alpha = p;
            alpha = 3 * p * k / len;
        }
        alpha = (alpha < 1.0) ? alpha : 1.0;
    }

}
*/

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

    #if VERBOSE > 1
        fprintf(stderr, ">>>> %d %lf %lf\n", up - low, gh1, gh);
        fflush(stderr);
    #endif

    return (gh <= gh1) ? gh : gh1;
}
