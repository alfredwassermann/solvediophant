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

#if defined(USEBLAS)
    #define BLAS 1
#else
    #define BLAS 0
#endif

#if BLAS
    #include "OpenBLASsub/common.h"
    #include "OpenBLASsub/cblas.h"
#endif

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

    #if 0
        swapvl = b[lattice->num_cols];
        for (i = lattice->num_cols; i > start; i--)
            b[i] = b[i - 1];
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

void insert_vector_long(lattice_t *lattice, long *u, int start, int end, int z) {
    long **b = lattice->basis_long;
    long *swap;
    int i, j, g;
    long q, ui;
    long hv;

    /* build new basis */
    for (j = 0; j < z; j++) {
        lattice->swap_long[j] = 0;
    }

    // Store new linear combination in lattice->swap
    for (i = start; i <= end; i++) {
        if (u[i] != 0) for (j = 0; j < z; j++) {
            lattice->swap_long[j] += b[i][j] * u[i];
        }
    }

    #if 0
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

void dual_insert_vector(lattice_t *lattice, long *u, int start, int end, int z, mpz_t hv) {
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

    #if 0
        swapvl = b[lattice->num_cols];
        for (i = lattice->num_cols; i > start; i--)
            b[i] = b[i - 1];
        b[start] = lattice->swap;
        lattice->swap = swapvl;
        lattice->num_cols++;
    #else
        g = start;
        while (u[g] == 0) g++;
        i = g + 1;
        while (labs(u[g]) > 1) {
            while (u[i] == 0) i++;
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

        // (b[g], b[g+1], ... , b[end]) -> (b[g+1], ... , b[end], b[g])
        swapvl = b[g];
        for (i = g; i < end; i++) {
            b[i] = b[i + 1];
        }
        b[end] = lattice->swap;
        coeffinit(b[end], z);

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

DOUBLE self_dual_bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, int p,
        int (*solutiontest)(lattice_t *lattice, int k)) {
    DOUBLE **R, *h_beta, *N;
    DOUBLE **H;
    DOUBLE r_tt;
    DOUBLE new_cj;
    DOUBLE lD;

    static mpz_t hv;
    int zaehler = 0;
    int h, i, last, h_end;
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

    fprintf(stderr, "\n######### SELFDUAL BKZ ########\n");
    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) {
        u[i] = 0;
    }

    lllalloc(&R, &h_beta, &N, &H, s, z);

    while (1) {
        fprintf(stderr, "Start tour\n");
        lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest);

        fprintf(stderr, "Primal\n");
        //for (start_block = 0; start_block + beta - 1 <= last; ++start_block) {
        //    end_block = start_block + beta - 1;
        for (start_block = 0; start_block < last; ++start_block) {
            end_block = start_block + beta - 1;
            if (end_block > last) end_block = last;

            new_cj = enumerate(lattice, R, u, s, start_block, end_block, delta, p);
            h = (start_block - 1 < 0) ? 0 : start_block - 1;
            h_end = (end_block + 1 <= last) ? end_block + 1 : last;

            r_tt = R[start_block][start_block];
            r_tt *= r_tt;
            if (delta * r_tt > new_cj) {
                fprintf(stderr, "primal enumerate successful %d %lf improvement: %lf\n",
                    start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
                fflush(stderr);
                insert_vector(lattice, u, start_block, end_block, z, hv);
                lllH(lattice, R, h_beta, H, h, 0, h_end, z, delta, CLASSIC_LLL, bit_size, solutiontest);
            } else {
                lllH(lattice, R, h_beta, H, h, h, h_end, z, 0.0, CLASSIC_LLL, bit_size, solutiontest);
            }
        }

        fprintf(stderr, "Dual\n");
        for (start_block = last - beta + 1; start_block > 0; --start_block) {
            end_block = start_block + beta - 1;

            new_cj = dual_enumerate(lattice, R, u, s, start_block, end_block, delta, p);
            h = (start_block - 1 < 0) ? 0 : start_block - 1;
            h_end = (end_block + 1 <= last) ? end_block + 1 : last;

            r_tt = 1.0 / R[end_block][end_block];
            r_tt *= r_tt;
            if (delta * r_tt > new_cj) {
                fprintf(stderr, "dual enumerate successful %d %lf improvement: %lf\n",
                    start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
                fflush(stderr);
                dual_insert_vector(lattice, u, start_block, end_block, z, hv);
                lllH(lattice, R, h_beta, H, h, 0, h_end, z, delta, CLASSIC_LLL, bit_size, solutiontest);
            } else {
                lllH(lattice, R, h_beta, H, h, h, h_end, z, 0.0, CLASSIC_LLL, bit_size, solutiontest);
            }
        }

        zaehler++;
        if (zaehler > 3) break;

    } /* end of |while| */

    lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest);

    lD = log_potential(R, s, z);

    #if 1
    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stderr);
    #endif

    lllfree(R, h_beta, N, H, s);
    free(u);
    mpz_clear(hv);

    return lD;
}

DOUBLE dual_bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, int p,
        int (*solutiontest)(lattice_t *lattice, int k)) {
    DOUBLE **R, *h_beta, *N;
    DOUBLE **H;
    DOUBLE r_tt;
    DOUBLE new_cj, new_cj2;
    DOUBLE lD;

    static mpz_t hv;
    int zaehler;
    int h, i, last, h_end;
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

    fprintf(stderr, "\n######### DUAL BKZ ########\n");
    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) {
        u[i] = 0;
    }

    lllalloc(&R, &h_beta, &N, &H, s, z);
    lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest);

    zaehler = -1;
    end_block = last + 1;
    //start_block = 0;
    while (zaehler < last) {
        end_block--;
        if (end_block < 2) {
            break;
        }
        h_end = (end_block < last) ? end_block + 1 : last;

        start_block = end_block - beta + 1;
        start_block = (start_block >= 0) ? start_block : 0;

        new_cj = dual_enumerate(lattice, R, u, s, start_block, end_block, delta, p);
        h = (start_block - 1 < 0) ? 0 : start_block - 1;

        r_tt = 1.0 / R[end_block][end_block];
        r_tt *= r_tt;
        if (delta * r_tt > new_cj) {
            fprintf(stderr, "dual enumerate successful %d %lf improvement: %lf\n",
                start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
            fflush(stderr);

            /* successful enumeration */
            dual_insert_vector(lattice, u, start_block, end_block, z, hv);
            //i = householder_column(lattice->basis, R, H, h_beta, end_block, end_block + 1, z, bit_size);
            lllH(lattice, R, h_beta, H, h, h, h_end, z, 0.0, CLASSIC_LLL, bit_size, solutiontest);
            i = end_block;

            new_cj2 = 1.0 / (R[i][i] * R[i][i]);
            if (FALSE && fabs(new_cj2 - new_cj) > EPSILON) {
                fprintf(stderr, "???????????????? We have a problem at %d: %lf %lf\n", i, new_cj2, new_cj);
                fflush(stderr);
                exit(1);
            }
            lllH(lattice, R, h_beta, H, h, 0, h_end, z, delta, CLASSIC_LLL, bit_size, solutiontest);
            //zaehler = -1;
            zaehler++;
        } else {
            zaehler++;
        }
    } /* end of |while| */

    lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest);

    lD = log_potential(R, s, z);

    #if 0
    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stderr);
    #endif

    lllfree(R, h_beta, N, H, s);
    free(u);
    mpz_clear(hv);

    return lD;
}

/**
 * Blockwise Korkine Zolotareff reduction
 */
DOUBLE bkz(lattice_t *lattice, int s, int z, DOUBLE delta, int beta, int p,
            int (*solutiontest)(lattice_t *lattice, int k), int (*solutiontest_long)(lattice_t *lattice, int k)) {
    DOUBLE **R, *h_beta, *N;
    DOUBLE **H;
    DOUBLE r_tt;
    DOUBLE new_cj, new_cj2;
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

    fprintf(stderr, "\n######### BKZ ########\n");
    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) {
        u[i] = 0;
    }

    lllalloc(&R, &h_beta, &N, &H, s, z);
    if (bit_size < 32) {
        copy_lattice_to_long(lattice);
        lllH_long(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest_long);
    } else {
        lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size, solutiontest);
    }

    start_block = zaehler = -1;
    //start_block = 0;
    while (zaehler < last) {
        start_block++;
        if (start_block == last) {
            start_block = 0;
        }
        if (start_block == 0) {
            //lllH(lattice, R, h_beta, H, 0, 0, s, z, delta, POT_LLL, bit_size);
        }

        end_block = start_block + beta - 1;
        end_block = (end_block < last) ? end_block : last;

        new_cj = enumerate(lattice, R, u, s, start_block, end_block, delta, p);
        h = (end_block + 1 < last) ? end_block + 1 : last;

        r_tt = R[start_block][start_block];
        r_tt *= r_tt;
        if (delta * r_tt > new_cj) {
            fprintf(stderr, "enumerate successful %d %lf improvement: %lf\n",
                start_block,  delta * r_tt - new_cj, new_cj / (delta * r_tt));
            fflush(stderr);

            /* successful enumeration */
            if (bit_size < 32) {
                insert_vector_long(lattice, u, start_block, end_block, z);
            } else {
                insert_vector(lattice, u, start_block, end_block, z, hv);
            }
            i = householder_column_long(lattice->basis_long, R, H, h_beta, start_block, start_block + 1, z, bit_size);
            new_cj2 = R[i][i] * R[i][i];
            if (fabs(new_cj2 - new_cj) > EPSILON) {
                fprintf(stderr, "???????????????? We have a problem: %lf %lf\n", new_cj2, new_cj);
                fflush(stdout);
            }

            if (bit_size < 32) {
                lllH_long(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, delta, CLASSIC_LLL, bit_size, solutiontest_long);
            } else {
                lllH(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, delta, CLASSIC_LLL, bit_size, solutiontest);
            }
            //lllH(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, 0.0, CLASSIC_LLL, bit_size, solutiontest);
            //lattice->num_cols--;

            //start_block = lllH(lattice, R, h_beta, H, start_block - 1, 0, h + 1, z, delta, CLASSIC_LLL, bit_size);
            //fprintf(stderr, "%d\n", start_block);

            zaehler = -1;
        } else {
            //fprintf(stderr, "enumerate: no improvement %d\n", zaehler);
            //fflush(stderr);
            if (h > 0) {
                if (bit_size < 32) {
                    lllH_long(lattice, R, h_beta, H, h - 1,  h - 1, h + 1, z, 0.0, CLASSIC_LLL, bit_size, solutiontest_long);
                } else {
                    lllH(lattice, R, h_beta, H, h - 1, h - 1, h + 1, z, 0.0, CLASSIC_LLL, bit_size, solutiontest);
                }
            }
            //start_block++;
            zaehler++;
        }
    } /* end of |while| */
    if (bit_size < 32) {
        copy_lattice_to_mpz(lattice);
    }

    lD = log_potential(R, s, z);

    #if 0
    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stderr);
    #endif

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

    //DOUBLE x;
    DOUBLE *y, *c;
    DOUBLE c_min;

    int i, j;
    int t, t_max;
    int found_improvement = 0;

    long *delta, *d, *v;
    DOUBLE *u_loc;
    int len, k;
    double alpha, radius;
    //DOUBLE *lambda_min;
    int SCHNITT = 20;

    c = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    y = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    delta = (long*)calloc(s+1,sizeof(long));
    d = (long*)calloc(s+1,sizeof(long));
    v = (long*)calloc(s+1,sizeof(long));
    u_loc = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    //lambda_min = (DOUBLE*)calloc(s+1,sizeof(DOUBLE));

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

    // Find minimum Eigen value
    /*
    i = start_block;
    lambda_min[i] = R[i][i] * R[i][i];
    for (i = start_block + 1; i <= end_block; ++i) {
        x = R[i][i] * R[i][i];
        lambda_min[i] = (x < lambda_min[i-1]) ? x : lambda_min[i-1];
    }
    */

    //t = t_max = end_block;
    for (t_max = start_block + 1; t_max <= end_block; t_max ++) {
        t = t_max;
        u_loc[t] = 1.0;

        while (t <= t_max) {
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
            // Use minimum Eigen value
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
                if (t > t_max) {
                    break;
                }
            }
            // next
            if (t == t_max) {
                //fprintf(stderr, "%d %ld %ld u=%lf v=%ld, %ld\n", t, d[t], delta[t], u_loc[t], v[t], v[t] + delta[t]);
                u_loc[t] += 1;
            } else {
                if (t < t_max) delta[t] *= -1.0;
                if (delta[t] * d[t] >= 0) delta[t] += d[t];
                u_loc[t] = v[t] + delta[t];
            }
        }
    }

    free(c);
    free(y);
    free(delta);
    free(d);
    free(v);
    free(u_loc);
    //free(lambda_min);

    if (!found_improvement) {
        c_min = R[start_block][start_block];
        c_min *= c_min;
    }
    return (c_min);
}


DOUBLE dual_enumerate(lattice_t *lattice, DOUBLE **R, long *u, int s,
                    int start_block, int end_block, DOUBLE improve_by, int p) {

    DOUBLE *y, *c, *a;
    DOUBLE c_min;

    int i, j;
    int t, t_min;
    int found_improvement = 0;

    long *delta, *d, *v;
    DOUBLE *u_loc;
    int len, k;
    double alpha, radius;
    int SCHNITT = 20;

    //fprintf(stderr, "-----------\n");
    c = (DOUBLE*)calloc(s+1, sizeof(DOUBLE));
    y = (DOUBLE*)calloc(s+1, sizeof(DOUBLE));
    a = (DOUBLE*)calloc(s+1, sizeof(DOUBLE));
    d = (long*)calloc(s+1, sizeof(long));
    v = (long*)calloc(s+1, sizeof(long));
    u_loc = (DOUBLE*)calloc(s+1, sizeof(DOUBLE));
    delta = (long*)calloc(s+1, sizeof(long));

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

        //fprintf(stderr, "LOOP %d %d %d\n", t_min, end_block, s+1);
        while (t >= t_min) {
            handle_signals(lattice, R);

            a[t] = u_loc[t] - y[t];
            /*
                c[t] = c[t - 1] + (a[t]/ R[t][t])^2
            */
            c[t] = a[t] / R[t][t];
            c[t] *= c[t];
            if (t > t_min) {
                c[t] += c[t - 1];
            }

            if (len <= SCHNITT) {
                alpha = 1.0;
            } else {
                k = t - start_block + 1;
                if (k > 2 * len / 4) {
                    alpha = 1.0;
                } else {
                    alpha = 0.5;
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

    free(delta);
    free(u_loc);
    free(v);
    free(d);
    free(a);
    free(y);
    free(c);

    if (!found_improvement) {
        c_min = 1.0 / R[end_block][end_block];
        c_min *= c_min;
    }
    return (c_min);
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

    #if VERBOSE > 1
        fprintf(stderr, ">>>> %d %lf %lf\n", up - low, gh1, gh);
        fflush(stderr);
    #endif

    return (gh <= gh1) ? gh : gh1;
}
