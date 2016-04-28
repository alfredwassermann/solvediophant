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
#include "OpenBLAS/common.h"
#include "OpenBLAS/cblas.h"

#define BLAS 1
#define DEEPINSERT 1
#define DEEPINSERT_CONST 100
#define VERBOSE 1

#define GIVENS 1
#define LASTLINESFACTOR "1000000" /* "100000000" */
#define EPSILON 0.000001      /* 0.0001  */
#define LLLCONST_LOW  0.7 /* 0.75*/
#define LLLCONST_HIGH 0.90    /* 0.99 */
#define LLLCONST_HIGHER 0.999
#define ETACONST 0.51

/**
 * global variables
 */
mpz_t matrix_factor;
mpz_t max_norm;
mpz_t max_norm_initial;
mpz_t max_up;
mpz_t dummy;

long nom, denom;
mpz_t lastlines_factor;

int system_rows, system_columns;
int lattice_rows, lattice_columns;
coeff_t **lattice;
int free_RHS;
int iszeroone;
mpz_t *upperbounds;
mpz_t upperbounds_max;
mpz_t upfac;

int *original_columns;
int no_original_columns;
int cut_after_coeff;
long stop_after_solutions;
long stop_after_loops;
long nosolutions;
//int iterate;
//int no_iterates;
//int bkz_beta, bkz_p;
int SILENT;
int nboundvars;

mpz_t soltest_u;
mpz_t soltest_s;
mpz_t soltest_upfac;

static FILE* fp;

/**
 * Inline functions
 */
#define ROUND(r) ceil(r-0.5)
#define put_to(i,j,val) mpz_set(lattice[i][j+1].c,val)
#define smult_lattice(i,j,factor) mpz_mul(lattice[i][j+1].c,lattice[i][j+1].c,factor)
#define get_entry(i,j) lattice[i][j+1].c

long diophant(gls_t *GLS, lll_params_t *LLL_params,
    mpz_t factor_input, mpz_t norm_input,
    int silent,
    long stop_after_sol_input,
    long stop_after_loops_input,
    int free_RHS_input,
    int cut_after, FILE* solfile) {

    int i,j;
    DOUBLE lD, lDnew;
    coeff_t *swap_vec;

    /**
     * Initialize some globals
     */
    mpz_init_set(matrix_factor,factor_input);
    mpz_init_set(max_norm,norm_input);
    mpz_init(lastlines_factor);
    mpz_init(upfac);

    mpz_init(soltest_u);
    mpz_init(soltest_s);
    mpz_init_set_ui(soltest_upfac, 1);

    /*
    if (LLL_params->iterate) {
        no_iterates = LLL_params->iterate_no;
    } else {
        bkz_beta = LLL_params->bkz.beta;
        bkz_p = LLL_params->bkz.p;
    }
    */
    SILENT = silent;
    stop_after_solutions = stop_after_sol_input;
    stop_after_loops = stop_after_loops_input;
    free_RHS = free_RHS_input;
    nom = 1;
    denom = 2;

    system_rows = GLS->num_rows;
    system_columns = GLS->num_cols;
    nboundvars = GLS->num_boundedvars;

#if BLAS
    //openblas_set_num_threads(8);
#endif

    /* In case, a time limit as been set by -time
       the execution is stopped and
       the number of solutions is printed
    */

    /**
     * set the lattice dimension;
     */
    lattice_rows = system_rows + system_columns + 1;
    lattice_columns = system_columns + 2;

    if (free_RHS) {
        lattice_rows++;
        lattice_columns++;
    } else {
        fprintf(stderr,"The RHS is fixed !\n");
        fflush(stderr);
    }
    cut_after_coeff = cut_after;

    /**
     * allocate memory
     */
    lattice = (coeff_t**)calloc(lattice_columns,sizeof(coeff_t*));
    for (j = 0; j < lattice_columns; j++) {
        lattice[j] = (coeff_t*)calloc(lattice_rows + 1, sizeof(coeff_t));
        for (i = 0; i <= lattice_rows; i++)
            mpz_init(lattice[j][i].c);
    }

    /**
     * read the system
     */
    for (j = 0; j < system_rows; j++) {
        for (i = 0; i < system_columns; i++) {
            mpz_mul(lattice[i][j+1].c, GLS->matrix[j][i], matrix_factor);
        }
        mpz_mul(lattice[system_columns][j+1].c, GLS->rhs[j], matrix_factor);
    }

    /**
     * handle upper bounds
     */
    mpz_init_set_si(upperbounds_max,1);
    iszeroone = 1;
    if (GLS->upperbounds == NULL) {
        fprintf(stderr, "No upper bounds: 0/1 variables are assumed \n"); fflush(stderr);
    } else {
        upperbounds = (mpz_t*)calloc(system_columns, sizeof(mpz_t));
        for (i = 0; i < system_columns; i++)
            mpz_init_set_si(upperbounds[i], 1);
        for (i = 0; i < nboundvars/*|system_columns|*/; i++) {
            mpz_set(upperbounds[i], GLS->upperbounds[i]);
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
    if (GLS->original_cols != NULL)
        no_original_columns = GLS->num_original_cols;
    else
        no_original_columns = GLS->num_cols;

    original_columns = (int*)calloc(GLS->num_original_cols, sizeof(int));

    if (GLS->original_cols != NULL)
        for (i = 0; i < no_original_columns; i++)
            original_columns[i] = GLS->original_cols[i];
    else {
        for (i = 0; i < no_original_columns; i++)
            original_columns[i] = 1;
        fprintf(stderr, "No preselected columns \n");
        fflush(stderr);
    }

    /**
     * append the other parts of lattice
     */
    for (j=system_rows;j<lattice_rows;j++) {
        mpz_mul_si(lattice[j-system_rows][j+1].c,max_norm, denom);
        mpz_mul_si(lattice[lattice_columns-2][j+1].c,max_norm, nom);
    }
    mpz_set(lattice[system_columns+free_RHS][lattice_rows].c, max_norm);

    if (free_RHS) {
        mpz_set_si(lattice[system_columns][lattice_rows-1].c, 1);
        mpz_set_si(lattice[system_columns+1][lattice_rows-1].c, 0);
    }
    mpz_set(lattice[system_columns+free_RHS][lattice_rows].c, max_norm);
    for (i=0;i<lattice_columns-1;i++) coeffinit(lattice[i], lattice_rows);

    /**
     * open solution file
     */
    fp = solfile;
    if (SILENT) fprintf(fp,"SILENT\n");
    fflush(fp);


#if 0
    printf("Before scaling\n");
    print_lattice();
#endif
    /**
     * scale lattice
     */
    mpz_init_set(max_norm_initial,max_norm);
    mpz_init_set_si(max_up,1);
    if (!iszeroone){
        for (j=0;j<nboundvars/*|system_columns|*/;j++) {
            if (mpz_sgn(upperbounds[j])!=0) {
                mpz_divexact(upfac,upperbounds_max,upperbounds[j]);
            } else {
                mpz_mul(upfac,upperbounds_max,upperbounds_max);
                mpz_mul_si(upfac,upfac,10000);
            }
            smult_lattice(j,j+system_rows, upfac );
            smult_lattice(system_columns+free_RHS,j+system_rows,upperbounds_max);
        }
        mpz_set(max_up,upperbounds_max);
        mpz_mul(max_norm,max_norm,max_up);
        if (free_RHS)
            smult_lattice(system_columns,lattice_rows-2,max_up);

        smult_lattice(system_columns+free_RHS,lattice_rows-1,max_up);
    }

#if 0
    printf("After scaling\n");
    print_lattice();
#endif

#if 1 // Do reduction
#if 0
    print_NTL_lattice();   /* Version for the NTL output */
    return 0;
#endif

    /**
     * permute lattice columns
     */
    swap_vec = lattice[lattice_columns-2];
    for (i = lattice_columns - 2; i > 0; i--)
        lattice[i] = lattice[i-1];
    lattice[0] = swap_vec;

#if 1
    printf("After permute\n");
    print_lattice();
#endif
    //shufflelattice();
    /**
     * first reduction
     */
    mpz_set_ui(lastlines_factor, 1);
    fprintf(stderr, "\n"); fflush(stderr);
    lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_LOW);

#if 1
    printf("After first reduction\n");
    print_lattice();
#endif
    /**
     * cut the lattice
     */
    if (cutlattice()) {
        fprintf(stderr, "First reduction successful\n"); fflush(stderr);
    } else {
        fprintf(stderr, "First reduction not successful\n"); fflush(stderr);
        return 0;
    }

#if 1
    printf("After cutting\n");
    print_lattice();
#endif

#if 1
    //shufflelattice();
    /**
     * second reduction
     */
    mpz_set_ui(lastlines_factor, 1);
    lll(lattice, lattice_columns-1, lattice_rows, LLLCONST_HIGH);
    fprintf(stderr, "Second reduction successful\n"); fflush(stderr);
#endif

#if 0
    printf("After second reduction\n");
    print_lattice();
#endif

#if 1  // Third reduction
    /**
     * scale last rows
     */
    mpz_set(lastlines_factor, LLL_params->scalelastlinefactor);
    for (i=0;i<lattice_columns;i++)
        mpz_mul(lattice[i][lattice_rows].c,lattice[i][lattice_rows].c,lastlines_factor);
    if (free_RHS)
        for (i=0;i<lattice_columns;i++)
            mpz_mul(lattice[i][lattice_rows-1].c,lattice[i][lattice_rows-1].c,lastlines_factor);

#if 0
    for (i=0;i<lattice_columns;i++) {
        for (j=0;j<40;j++)
            mpz_mul_ui(lattice[i][j+1].c,lattice[i][j+1].c, 9);
    }
#endif

    /**
     * third reduction
     */
    fprintf(stderr, "\n"); fflush(stderr);
    if (LLL_params->iterate) {
        iteratedlll(lattice, lattice_columns-1, lattice_rows, LLL_params->iterate_no, LLLCONST_HIGH);
    } else {
        shufflelattice();
        lDnew = bkz(lattice, lattice_columns, lattice_rows, LLLCONST_HIGHER,
                        LLL_params->bkz.beta, LLL_params->bkz.p);

        i = 0;
        do {
            lD = lDnew;
            shufflelattice();
            lDnew = bkz(lattice,lattice_columns,lattice_rows,LLLCONST_HIGH,
                        LLL_params->bkz.beta, LLL_params->bkz.p);
            printf("%0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD - lDnew);
            i++;
        }
        while (i < 1 && fabs(lDnew - lD) > 0.01);
    }
    fprintf(stderr, "Third reduction successful\n"); fflush(stderr);

    /* undo scaling of last rows */
    for (i = 0; i < lattice_columns; i++)
        mpz_divexact(lattice[i][lattice_rows].c, lattice[i][lattice_rows].c, lastlines_factor);
    if (free_RHS)
        for (i = 0; i < lattice_columns; i++)
            mpz_divexact(lattice[i][lattice_rows-1].c,
                lattice[i][lattice_rows-1].c,
                lastlines_factor);
#endif // Third reduction
#else
    read_NTL_lattice();
#endif // Do reduction

#if 1
    printf("Before enumeration\n");
    /*|print_NTL_lattice();|*/   /* Version for the NTL output */
    print_lattice();
#endif

    /**
     * explicit enumeration
     */
    fprintf(stderr, "\n"); fflush(stderr);
    nosolutions = explicit_enumeration(lattice,lattice_columns-1,lattice_rows);

    /**
     * close solution file;
     */
    print_num_solutions(nosolutions);

    /**
     * free multiprecision memory
     */
    mpz_clear(matrix_factor);
    mpz_clear(max_norm);
    mpz_clear(lastlines_factor);
    mpz_clear(upfac);
    mpz_clear(max_norm_initial);
    mpz_clear(max_up);
    mpz_clear(soltest_u);
    mpz_clear(soltest_s);
    mpz_clear(soltest_upfac);
    mpz_clear(upperbounds_max);

    for(j=0;j<lattice_columns;j++) {
        for (i=0;i<=lattice_rows;i++) mpz_clear(lattice[j][i].c);
    }
    free(lattice);
    if (upperbounds!=NULL) {
        for (i=0;i<system_columns;i++) mpz_clear(upperbounds[i]);
        free(upperbounds);
    }

    return nosolutions;
}

/**
 * Basic subroutines
 */
 void print_num_solutions(long num_solutions) {
     if (SILENT) fprintf(fp, "%ld solutions\n", nosolutions);
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
void print_lattice() {
    int i, j;
    for (i = 0; i < lattice_columns; i++) {
        for (j = 0; j < lattice_rows; j++) {
            mpz_out_str(NULL,10,get_entry(i,j));
            printf(" ");
        }
        printf("\n");
    }
    printf("\n");
    fflush(stdout);
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

int cutlattice() {
    int j, i, flag;

    /**
     * delete unnecessary columns
     */
    j=0;
    do {
        if (lattice[j][0].p > system_rows)
            j++;
        else {
            for (i = j + 1; i < lattice_columns; i++)
                lattice[i-1] = lattice[i];
            lattice_columns--;
        }
    } while (j < lattice_columns - 1);


    /**
     * test for right hand side columns
     */
    flag = 0;
    for (i = 0; i < lattice_columns; i++) if (mpz_sgn(get_entry(i,lattice_rows-1)) != 0) {
        flag = 1;
        break;
    }
    if (flag == 0) {
        fprintf(stderr, "Nonhomogenous solution not possible.\n"); fflush(stderr);
        exit(2);

        return 0;  /* Just for the compiler */
    }

    /* Now the rows are deleted. */
    for (j = 0; j < lattice_columns; j++)  {
       if (nboundvars == 0) {
            for (i = system_rows; i < lattice_rows; i++)
                put_to(j,i-system_rows,get_entry(j,i));
        } else {
            for (i = system_rows; i < system_rows + nboundvars; i++)
                put_to(j,i-system_rows,get_entry(j,i));
            for (i = system_rows+system_columns; i < lattice_rows; i++)
                put_to(j,i-system_rows-system_columns+nboundvars,get_entry(j,i));
        }
    }
    lattice_rows -= system_rows;
    lattice_rows -= (system_columns-nboundvars);

    for (j = 0; j < lattice_columns; j++) coeffinit(lattice[j],lattice_rows);

    return 1;
}

int solutiontest(int position) {
    int i,j;
    int low, up;
    int end;

    /* test the last two rows */
    if (mpz_cmpabs(get_entry(position,lattice_rows-1),max_norm)!=0) return 0;
    if (mpz_sgn(get_entry(position,lattice_rows-1-free_RHS))==0) return 0;

    /* test, if column is a solution */
    low = 0;
    up = lattice_rows-1-free_RHS;
    if (lattice_columns == system_columns+2+free_RHS) {
        for (i=0;i<system_rows;i++)
            if (mpz_sgn(get_entry(position,i))!=0) return 0;
        low = system_rows;
    }

    if (iszeroone) {
        for (i=low;i<up;i++) {
            if (mpz_cmpabs(get_entry(position,i),max_norm)!=0) return 0;
        }
    } else {
        for (i=low;i<up;i++) {
            if (mpz_cmpabs(get_entry(position,i),max_norm)>0) return 0;
        }
    }


    mpz_set_si(upfac,1);
    mpz_divexact(soltest_s,get_entry(position,lattice_rows-1),lastlines_factor);

    /* write a solution with blanks */
    i = low;
    if (cut_after_coeff == -1) {
        end=no_original_columns;
    } else {
        end=cut_after_coeff;
    }

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
            mpz_set(soltest_u,get_entry(position,i));
            mpz_sub(soltest_u,soltest_u,soltest_s);
            mpz_divexact(soltest_u,soltest_u,max_norm_initial);
            mpz_divexact(soltest_u,soltest_u,soltest_upfac);
            mpz_divexact_ui(soltest_u,soltest_u,denom);
            mpz_abs(soltest_u,soltest_u);
            i++;
        }
        mpz_out_str(NULL,10,soltest_u);
        printf(" ");
        if (stop_after_solutions==1) {
            mpz_out_str(fp,10,soltest_u);
            fprintf(fp," ");
        }
    }
    if (free_RHS) {
        mpz_divexact(soltest_u,get_entry(position,up),max_up);
        mpz_divexact(soltest_u,soltest_u,lastlines_factor);
        mpz_abs(soltest_u,soltest_u);
        printf(" L = ");
        mpz_out_str(NULL,10,soltest_u);
    }
    printf("\n");
    fflush(stdout);

    /* test if one solution is enough */
    if (stop_after_solutions == 1) {
        fprintf(fp,"\n");

        fprintf(stderr, "Stopped in phase 1 after finding a random solution\n");
        exit(8);
    }

    return 1;
}

/**
 *  Lattice basis reduction algorithms
 */

#define TWOTAUHALF 67108864.0 /* $2^{\tau/2}$*/

int lllfp(coeff_t **b, DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs,
            int start, int s, int z, DOUBLE delta) {

    int i, j, k;
    DOUBLE ss;

    int Fc, Fr;
    DOUBLE mus, cc;
    mpz_t musvl;
    mpz_t hv;
    DOUBLE *swapd;

    coeff_t *swapvl;

#if VERBOSE > 3
    int counter;
#endif

    int Fi;

    mpz_init(musvl);
    mpz_init(hv);

    /* Test for trivial cases. */
    if ((z <= 1) || (s <= 1)) {
        fprintf(stderr, "Wrong dimensions in lllfp\n");
        fflush(stderr);
        return(0);
    }

    k = (start > 1) ? start : 1;

    /* first step: compute norms
    Step 1:
    We begin with stage $k=1$. Schnorr's original algorithm runs
    from $k=2$ to $k=s$. Here we have $k=1,\ldots,s-1$.

    Initially the $b_i$'s are rounded and then their norms,
    $N_i = ||b'_i||$, are computed.

    Then we do steps 2 -- 4 while  $k\leq s$.
    */
    if (k < 1) k = 1;
    for (i = k - 1; i < s; ++i) {
        ss = 0.0;
        for (j = 0; j < z; ++j) {
            bs[i][j] = (DOUBLE)mpz_get_d(b[i][j+1].c);
            ss += bs[i][j]*bs[i][j];
        }
        N[i] = SQRT(ss);
    }

#if VERBOSE > 3
    counter = 0;
#endif

    /* The main loop */
    while (k < s) {
#if VERBOSE > 3
        if ((counter % 500) == 0) {
            fprintf(stderr, "LLL: %d k:%d\n", counter, k);
            fflush(stderr);
        }
        counter++;
#endif

        /* second step: orthogonalization of $b_k$
          Step 2:
          computation of $\mu_{k,1},\ldots,\mu_{k,k-1}$ and $c_k = ||\hat{b}_k||^2.$
          This is done with Gram Schmidt Orthogonalization.
        */
        if (k == 1) {
            c[0] = N[0] * N[0];
        }

        c[k] = N[k] * N[k];
        for (j = 0; j < k; j++) {
            ss = scalarproductfp(bs[k], bs[j], z);
            if (fabs(ss) < N[k] * N[j] / TWOTAUHALF) {
                ss = (DOUBLE)scalarproductlfp(b[k],b[j]);
            }
            for (i = 0; i < j; i++) {
                ss -= mu[j][i] * mu[k][i] * c[i];
            }
            mu[k][j] = ss / c[j];
            c[k] -= ss * mu[k][j];
        }

if (c[k] < EPSILON) {
    fprintf(stderr, "c[%d] is very small: %lf\n", k, c[k]);
}

        /* third step: size reduction of $b_k$ */
        Fc = Fr = 0;
        for (j = k - 1; j >= 0; j--) {
            if (fabs(mu[k][j]) > ETACONST) {
                /* round the Gram Schmidt coefficient */
                mus = ROUND(mu[k][j]);
                mpz_set_d(musvl, mus);
                if (fabs(mus) > TWOTAUHALF) {
                    Fc = 1;
                }

                Fr = 1;
                /* set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$ */
                size_reduction(b, mu, musvl, mus, k, j);

                if (Fc == 1)
                    fprintf(stderr, "Problems in rounding mu, step back to k=%d; %lf %lf\n", k-1, mus, mu[k][j]);

                mu[k][j] -= mus;
                solutiontest(k);
            }
        }

        /* recompute $N_k$ */
        N[k] = 0.0;
        for(i = 0; i < z; i++) {
            bs[k][i] = (DOUBLE)mpz_get_d(b[k][i+1].c);
            N[k] += bs[k][i]*bs[k][i];
        }
        N[k] = SQRT(N[k]);

        /* Step back because of floating point instability */
        if (Fc == 1) {
            k = (k - 1 > 1) ? k - 1 : 1;    /* $k=\max(k-1,1)$ */
        } else {
            /* test for linear dependencies */

            if (N[k] <- EPSILON) {
                fprintf(stderr, "Nk negativ! contact the author.\n");
                fflush(stderr);
                exit(1);
            }
            /*
            Before going to step 4 we test if $b_k$ is linear dependent.
            This is the case if $||N_k||<{1\over 2}$.
            If we found a linear dependent vector $b_k$,
            we shift $b_k$ to the last column of the
            matrix and restart |lllfp| with $s:= s-1$.
            */
            if (N[k] < 0.5) {
               swapvl = b[k];
               ss = N[k];
               swapd = bs[k];
               for (i = k + 1; i < s; i++) {
                    b[i-1] = b[i];
                    N[i-1] = N[i];
                    bs[i-1] = bs[i];
               }
               b[s-1] = swapvl;
               N[s-1] = ss;
               bs[s-1] = swapd;
               s = s-1;
               k = 1;
               continue;
            }

        }

#if defined(DEEPINSERT)
        /* fourth step: deepinsert columns */
        cc = N[k] * N[k];
        j = 0;
        Fi = 0;
        while (j < k) {
#if 1
            if ((j > DEEPINSERT_CONST && j < k - DEEPINSERT_CONST) || delta * c[j] <= cc) {
#else
            if (delta * c[j] <= cc) {
#endif
                cc -= mu[k][j]*mu[k][j]*c[j];
                j++;
            } else {
                swapvl = b[k];
                ss = N[k];
                swapd = bs[k];
                for(i=k-1;i>=j;i--) {
                    b[i+1] = b[i];
                    N[i+1] = N[i];
                    bs[i+1] = bs[i];
                }
                b[j] = swapvl;
                N[j] = ss;
                bs[j] = swapd;

                Fi = 1;
                break;
            }
        }
        if (Fi == 1)
          k = (j-1 > 1) ? j-1 : 1;       /* $k = \max(j-1,1)$ */
        else {
            k++;
        }

#else
        /* fourth step: swap columns */
        if (delta * c[k-1] > c[k] + mu[k][k-1] * mu[k][k-1] * c[k-1]) {
            swapvl = b[k];
            b[k] = b[k-1];
            b[k-1] = swapvl;
            ss = N[k];
            N[k] = N[k-1];
            N[k-1] = ss;
            swapd = bs[k];
            bs[k] = bs[k-1];
            bs[k-1] = swapd;

            k = (k-1 > 1) ? k-1 : 1; /* $k = \max(k-1,1)$ */
        } else
            k++;

#endif
    }
    mpz_clear(hv);
    mpz_clear(musvl);
    return(1);

}

int lllHfp(coeff_t **b, DOUBLE **R, DOUBLE *c, DOUBLE *N, DOUBLE **H,
            int start, int s, int z, DOUBLE delta) {

    int i, j, k;
    DOUBLE ss, x;
    DOUBLE zeta;
    DOUBLE rkk, rii;
    DOUBLE beta[32768];
    DOUBLE w;
    DOUBLE bb;
    DOUBLE norm;

    DOUBLE mus;
    mpz_t musvl;
    mpz_t sum_mu;
    mpz_t hv;

    coeff_t *swapvl;

#if VERBOSE > 3
    int counter;
#endif

fprintf(stderr, "------------------------NEW LLLHfp-----------------------------\n");
fflush(stderr);

    mpz_init(musvl);
    mpz_init(hv);
    mpz_init(sum_mu);

    /* Test for trivial cases. */
    if ((z <= 1) || (s <= 1)) {
        fprintf(stderr, "Wrong dimensions in lllfp\n");
        fflush(stderr);
        return(0);
    }

    k = (start > 0) ? start : 0;
    //if (k < 1) k = 1;
    k = 0;

#if VERBOSE > 3
    counter = 0;
#endif

    /* The main loop */
    while (k < s) {
#if VERBOSE > 3
        if ((counter % 500) == 0) {
            fprintf(stderr, "LLL_H: %d k:%d\n", counter, k);
            fflush(stderr);
        }
        counter++;
#endif

//fprintf(stderr, "\nk %d\n", k);
start_tricol:

        /* Recompute column k of R */
        for (j = 0; j < z; ++j) {
            R[k][j] = (DOUBLE)mpz_get_d(b[k][j+1].c);
        }

        for (i = 0; i < k; ++i) {
            for (j = i, w = 0.0; j < z; ++j) {
                w += R[k][j] * H[i][j];
            }
            for (j = 0; j < z; ++j) {
                R[k][j] -= w * beta[i] * H[i][j];
            }
        }
        /* Now, R[k] is updated. */

        // Norm of pivot column
        for (j = k, norm = 0.0; j < z; ++j) {
            norm += R[k][j] * R[k][j];
        }
        norm = sqrt(norm);

        for (j = k; j < z; ++j) {
            H[k][j] = R[k][j] / norm;
        }
        H[k][k] += (R[k][k] >= 0.0) ? 1 : -1;

        bb = 1.0 / (1.0 + abs(R[k][k]) / norm);
        for (j = k, ss = 0.0; j < z; ++j) {
            ss += H[k][j] * H[k][j];
        }
        beta[k] = 2.0 / ss;

        for (j = k, w = 0.0; j < z; ++j) {
            w += R[k][j] * H[k][j];
        }
//fprintf(stderr, "beta %lf, beta_s %lf, w %lf\n", beta[k], bb, w);

        for (j = k; j < z; ++j) {
            R[k][j] -= beta[k] * w * H[k][j];
        }

#if 0
fprintf(stderr, "H\n");
for (i = 0; i <=k; i++) {
    for (j = 0; j < z; ++j) {
        fprintf(stderr, "%0.4lf ", H[i][j]);
    }
    fprintf(stderr, "\n");
}

fprintf(stderr, "R\n");
for (i = 0; i <=k; i++) {
    for (j = 0; j < z; ++j) {
        fprintf(stderr, "%0.4lf ", R[i][j]);
    }
    fprintf(stderr, "\n");
}
/*
for (i = 0; i < z; i++) {
    for (j = 0; j < z; ++j) {
        if (i == j)
            fprintf(stderr, "%0.4lf ", 1 - beta[k] * H[k][i] * H[k][j]);
        else
            fprintf(stderr, "%0.4lf ", 0 - beta[k] * H[k][i] * H[k][j]);
    }
    fprintf(stderr, "\n");
}
*/
#endif

//if (k > 1)
//    exit(1);
        /* third step: size reduction of $b_k$ */
        mpz_set_si(sum_mu, 0);
        for (j = k - 1; j >= 0; j--) {
            ss = R[k][j] / R[j][j];
            if (fabs(ss) > ETACONST) {
                mus = ROUND(ss);
//fprintf(stderr, "mu %lf\n", mus);
                mpz_set_d(musvl, mus);
                mpz_add_ui(sum_mu, sum_mu, (unsigned long)abs(mus));

                /* set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$ */
                size_reduction(b, R, musvl, mus, k, j);
                solutiontest(k);
            }
        }
        if (mpz_cmp_si(sum_mu, 0) != 0) {
//            fprintf(stderr, "REDO tricol\n");
            goto start_tricol;
        }


            /*
            Before going to step 4 we test if $b_k$ is linear dependent.
            This is the case if $||N_k||<{1\over 2}$.
            If we found a linear dependent vector $b_k$,
            we shift $b_k$ to the last column of the
            matrix and restart |lllfp| with $s:= s-1$.
            */

        /* fourth step: swap columns */
        if (k > 0 &&
            delta * R[k-1][k-1]*R[k-1][k-1] > R[k][k-1]*R[k][k-1] + R[k][k]*R[k][k]) {

fprintf(stderr, "SWAP %d %d\n", k-1, k);
            swapvl = b[k];
            b[k] = b[k-1];
            b[k-1] = swapvl;

            //k = (k-1 > 1) ? k-1 : 1; /* $k = \max(k-1,1)$ */
            k--;
        } else {
            k++;
        }
    }
    mpz_clear(hv);
    mpz_clear(musvl);
    return(1);

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
    return cblas_ddot(n,v,1,w,1);
#else
    DOUBLE r;
    int i;
    r = 0.0;
    for (i = n - 1; i >= 0; i--) r += v[i] * w[i];
    return r;
#endif
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
        for (i = 0; i < j; i++) mu[k][i] -= mu[j][i];
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
        for (i = 0; i < j; i++) mu[k][i] += mu[j][i];
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
#if 0
        daxpy(j,-mus,mu[k],1,mu[j],1);
#endif
        for (i = 0; i < j; i++) mu[k][i] -= mu[j][i] * mus;

    }
}

int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N,  DOUBLE ***bs, int s, int z) {
    int i, m;

    if ((z < 1) || (s < 1)) return 0;

    (*c) = (DOUBLE*)calloc(s, sizeof(DOUBLE));
    (*N) = (DOUBLE*)calloc(s, sizeof(DOUBLE));
    (*mu) = (DOUBLE**)calloc(z, sizeof(DOUBLE*));
    for (i = 0; i < s; i++)
        (*mu)[i] = (DOUBLE*)calloc(z, sizeof(DOUBLE));

    m = (z > s) ? z : s;
    (*bs) = (DOUBLE**)calloc(m,sizeof(DOUBLE*));

    for (i = 0; i < m; i++) (*bs)[i] = (DOUBLE*)calloc(z, sizeof(DOUBLE));

    return 1;
}

int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s) {
    int i;

    for (i = 0; i < s; ++i) free(bs[i]);
    free(bs);
    for (i = 0; i < s; ++i) free(mu[i]);
    free(mu);
    free(N);
    free(c);

    return 1;
}

double logD(coeff_t **lattice, DOUBLE *c, int s, int z) {
    double d = 0.0;
    int i;

    for (i = 0; i < s; i++) {
        d += log(c[i]) * (s - i);
    }
    d *= 0.5;
    return d;
}

double orthogonal_defect(coeff_t **lattice, DOUBLE *c, int s, int z) {
    double defect = 0.0;

#if 0
    int i;
    for (i=0;i<s;i++) defect += log((double)normfp(lattice[i]))
        - log((double)c[i]);
#endif

    defect /= 2.0;
    return defect;
}

/**
 * LLL variants
 */
void lll(coeff_t **b, int s, int z, DOUBLE quality) {
    DOUBLE **mu;
    DOUBLE *c;
    DOUBLE *N;
    DOUBLE **bs;
    int r;

    lllalloc(&mu, &c, &N, &bs, s, z);
    r = lllHfp(b, mu, c, N, bs, 1, s, z, quality);
    lllfree(mu, c, N, bs, s);

    return;
}

DOUBLE iteratedlll(coeff_t **b, int s, int z, int no_iterates, DOUBLE quality) {
    DOUBLE **mu;
    DOUBLE *c;
    DOUBLE *N;
    DOUBLE **bs;
    int r, l, i, j, runs;
    coeff_t *swapvl;
    DOUBLE lD;

    lllalloc(&mu,&c,&N,&bs,s,z);
    r = lllHfp(b, mu, c, N, bs, 1, s, z, quality);

    lD = logD(b,c,s,z);
    fprintf(stderr, "   log(D)= %f\n", lD);
    fflush(stderr);

    for (runs = 1; runs < no_iterates; runs++) {
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

        r = lllfp(b,mu,c,N,bs,1,s,z,quality);
        lD = logD(b,c,s,z);
        fprintf(stderr, "%d: log(D)= %f\n", runs, lD);
        fflush(stdout);
    }

    lllfree(mu, c, N, bs, s);

    return lD;
}

/**
 * Blockwise Korkine Zolotareff reduction
 */
DOUBLE bkz(coeff_t **b, int s, int z, DOUBLE delta, int beta, int p) {
    DOUBLE **mu, *c, *N;
    DOUBLE **bs;
    static mpz_t hv;
    int zaehler;
    int h,i,last;
    int start_block, end_block;
    long *u;
    DOUBLE new_cj;
    DOUBLE lD;

    int g, ui, q, j;
    coeff_t *swapvl;

    mpz_init(hv);

    last = s - 2;    /* |last| points to the last nonzero vector of the lattice.*/
    if (last < 1) {
        printf("BKZ: the number of basis vectors is too small.\n");
        printf("Probably the number of rows is less or equal");
        printf(" to number of columns in the original system\n");
        printf("Maybe you have to increase c0 (the first parameter)!\n");

        mpz_clear(hv);
        return 0.0;
    }

    u = (long*)calloc(s, sizeof(long));
    for (i = 0; i < s; i++) u[i] = 0;

    lllalloc(&mu,&c,&N,&bs,s,z);
    lllfp(b,mu,c,N,bs,1,s,z,delta);

    start_block = zaehler = -1;
    while (zaehler < last) {
        start_block++;
        if (start_block == last) start_block = 0;
        end_block = (start_block+beta-1 < last) ? start_block+beta-1 : last;
        new_cj = enumerate(mu,c,u,s,start_block,end_block,p);

        /* The exhaustive enumeration. */
        h = (end_block + 1 < last) ? end_block + 1 : last;

        if (delta * c[start_block] > new_cj) {
            /* successful enumeration */
            /* build new basis */
            for (j = 1; j <= z; j++) mpz_set_si(b[last+1][j].c, 0);
            for (i = start_block; i <= end_block; i++) {
                if (u[i] != 0) for(j = 1; j <= z; j++) {
                    if (u[i] > 0) {
                        mpz_addmul_ui(b[last+1][j].c,b[i][j].c,u[i]);
                    } else {
                        mpz_submul_ui(b[last+1][j].c,b[i][j].c,-u[i]);
                    }
                }
            }
            g = end_block;
            while (u[g] == 0) g--;

            i = g - 1;
            while (labs(u[g]) > 1) {
                while (u[i] == 0) i--;
                q = (int)ROUND((1.0 * u[g]) / u[i]);
                ui = u[i];
                u[i] = u[g] - q*u[i];
                u[g] = ui;

                for (j = 1; j <= z; j++) {
                    mpz_set(hv,b[g][j].c);
                    mpz_mul_si(b[g][j].c,b[g][j].c,(long)q);
                    mpz_add(b[g][j].c,b[g][j].c,b[i][j].c);
                    mpz_set(b[i][j].c,hv);
                }
                coeffinit(b[g],z);
                coeffinit(b[i],z);
            }

            swapvl = b[g];
            for (i = g; i > start_block; i--) b[i] = b[i-1];
            b[start_block] = b[last+1];
            coeffinit(b[start_block],z);

            b[last+1] = swapvl;
            for (j = 1; j <= z; j++) mpz_set_si(b[last+1][j].c,0);
            coeffinit(b[last+1],z);

            lllfp(b,mu,c,N,bs,start_block-1,h+1,z,delta);

            if (N[h]<-EPSILON) {
                fprintf(stderr,"NN negativ\n");
                fflush(stderr);
                printf("NN negativ\n");
                fflush(stdout);
                exit(1);
            }

            zaehler = -1;
        } else {
            if (h > 0) {
                lllfp(b,mu,c,N,bs,h-2,h+1,z,delta);   /* For some unkown reason we have to
                                                        use $h-2$ as |start|. */
            }
            zaehler++;
        }
    } /* end of |while| */

    lD = logD(b,c,s-1,z);

    fprintf(stderr, "bkz: log(D)= %f\n", lD);
    fflush(stdout);
    lllfree(mu,c,N,bs,s);
    free(u);
    mpz_clear(hv);

    return lD;
}

/**
 * Pruned Gauss-Enumeration.
 */
DOUBLE enumerate(DOUBLE **mu, DOUBLE *c, long *u, int s, int start_block, int end_block, int p) {
    DOUBLE cd, dum;
    DOUBLE *y, *cs, *eta;

    DOUBLE **sigma;
    int *r;

    long *us, *delta, *d, *v;
    int t,i,t_up, len;
    double alpha;
    int tmax;
    static DOUBLE pi = 3.141592653589793238462643383;
    static DOUBLE e = 2.718281828459045235360287471352662497757247093;
    int SCHNITT = 40;

    if (c[start_block]<=EPSILON) {
        fprintf(stderr, "Hier ist was faul! start_block=%d %f\n", start_block, (double)c[start_block]);
        fflush(stderr);
        exit(1);
    }

    us=(long*)calloc(s+1,sizeof(long));
    cs=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    y=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    delta=(long*)calloc(s+1,sizeof(long));
    d=(long*)calloc(s+1,sizeof(long));
    eta=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
    v=(long*)calloc(s+1,sizeof(long));

    sigma =(DOUBLE**)calloc(s,sizeof(DOUBLE*));
    r =(int*)calloc(s+1,sizeof(int));
    for (i = 0; i < s; i++) {
        sigma[i] =(DOUBLE*)calloc(s,sizeof(DOUBLE));
        r[i] = i-1;
    }

    len = end_block + 1 - start_block;
    for (i = start_block; i <= end_block + 1; i++) {
        cs[i] = y[i] = 0.0;
        u[i] = us[i] = v[i] = delta[i] = 0;
        d[i] = 1;
    }
    us[start_block] = u[start_block] = 1;
    cd = c[start_block];

    t = tmax = start_block;      /* Now we start from $t=|start_block|$ instead of
     $t=|end_block|$. */

    /* precompute $\eta$ */
    eta[start_block] = 0.0;
    if (end_block-start_block <= SCHNITT) {
        for (i = start_block + 1; i <= end_block; i++) eta[i] = 0.0;
    } else {
        dum = log(c[start_block]);
        /*
            Hoerners version of the Gaussian volume heuristics.
        */
        dum = log(c[start_block]);
        for (i = start_block + 1; i <= end_block; i++) {
            t_up = i - start_block;
            eta[i] = 0.5*t_up*exp( (log(pi*t_up)-2.0*p*log(2.0)+dum)/t_up )/(pi*e);
            if (i < end_block) dum += log(c[i]);
        }
    }


    while (t <= end_block) {
        /* the block search loop */
        dum = us[t] + y[t];
        cs[t] = cs[t+1] + dum*dum*c[t];

        if (len <= SCHNITT) {
            alpha = 1.0;
        } else {
            alpha = sqrt(1.20 * (end_block + 1 - t) / len);
            if (alpha >= 1.0) alpha = 1.0;
        }
        alpha *= cd;

        if (cs[t] < alpha - EPSILON) {
           if (t > start_block) {
               t--;
               if (r[t+1] > r[t]) r[t] = r[t+1];

               delta[t] = 0;
               for (i = r[t+1]; i > t; i--) sigma[i][t] = sigma[i+1][t] + us[i]*mu[i][t];

               y[t] = sigma[t+1][t]; /*|dum;|*/
               us[t] = v[t] = (long)(ROUND(-y[t]));
               d[t] = (v[t] > -y[t]) ? -1 : 1;
           } else {
               cd = cs[start_block];
               for (i = start_block; i <= end_block; i++) u[i] = us[i];
               goto nextstep;
           }
       } else {
           t++;
           r[t] = t;
nextstep:
           if (tmax < t) tmax = t;
           if (t < tmax) delta[t] = -delta[t];
           if (delta[t] * d[t] >= 0) delta[t] += d[t];
           us[t] = v[t] + delta[t];
       }

    }
    free (us);
    free (cs);
    free (y);
    free (delta);
    free (d);
    free (eta);
    free (v);
    for (i = s - 1; i >= 0; i--) {
       free(sigma[i]);
    }
    free(sigma);
    free(r);
    return (cd);
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
    DOUBLE **fipo_UB, **fipo_LB;
#endif
/*|mpz_t *upb,*lowb;|*/
long fipo_success;
DOUBLE dum1, dum2;

long only_zeros_no, only_zeros_success, hoelder_no, hoelder_success;
long cs_success;
long N_success;

DOUBLE explicit_enumeration(coeff_t **lattice, int columns, int rows) {
    /* local variables for |explicit_enumeration() */
    /*|__attribute((aligned(16)))|*/

    int level,level_max;
    int i,j,l;
    long loops;

    DOUBLE *y, *cs, *us;

    long *delta, *d, *eta;
    long *v;
    int *first_nonzero, *first_nonzero_in_column, *firstp;

    DOUBLE *N, **mu, *c, **w, **bd, **mu_trans;

    DOUBLE Fd, Fq, Fqeps;
    DOUBLE *dum;
    DOUBLE tmp;
    // coeff_t *swap_vec;

    int isSideStep = 0;
    DOUBLE stepWidth = 0.0;
    DOUBLE olddum;

#if defined(FINCKEPOHST)
    DOUBLE *fipo;
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
    fipo=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    muinv=(DOUBLE**)calloc(columns,sizeof(DOUBLE*));
    for(i = 0; i < columns; ++i)
        muinv[i] = (DOUBLE*)calloc(rows,sizeof(DOUBLE));

    fipo_LB=(DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
    fipo_UB=(DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
    for (i = 0; i <= columns; ++i) {
        fipo_LB[i] = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
        fipo_UB[i] = (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    }
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

    /* count nonzero entries in the last rows(s) */
    if (free_RHS) {
        i=0;
        for (j = columns - 1; j >= 0; j--) if (mpz_sgn(get_entry(j,rows-2)) != 0)
            i++;
        fprintf(stderr, "Number of nonzero entries in the second last row: %d\n", i);
        fflush(stderr);
    }

    i = 0;
    for (j = columns - 1; j >= 0; j--) if (mpz_sgn(get_entry(j,rows-1)) !=0 )
        i++;
    fprintf(stderr, "Number of nonzero entries in the last row: %d\n", i);
    fflush(stderr);

#if 0
    @<sort lattice columns@>;
#endif

    /* set the simple pruning bounds */
    Fq = (DOUBLE)mpz_get_d(max_norm);
    Fd = (rows*Fq*Fq)*(1.0+EPSILON);
    Fqeps = (1.0 + EPSILON) * Fq;        /* Used in prune() */
#if VERBOSE > 0
    fprintf(stderr, "Fq: %f\n", (double)Fq);
    fprintf(stderr, "Fd: %f\n", (double)Fd);
    fflush(stderr);
#endif

    /* orthogonalize the basis */
#if GIVENS
    givens(lattice,columns,rows,mu,bd,c,N,Fq);
#else
    gramschmidt(lattice,columns,rows,mu,bd,c,N,Fq);
#endif

    /* compute $mu^\top$, the transpose of $mu$. */
    for (i = 0; i < columns; i++)
        for (j = 0; j < columns; j++)
            mu_trans[j][i] = mu[i][j];

#if 0
    basis2poly();
#endif

#if FINCKEPOHST
    /* determine Fincke-Pohst bounds */
    fipo_success = 0;
    inverse(mu, muinv, columns);

#if VERBOSE > -1
    fprintf(stderr, "\nFincke-Pohst bounds:\n");
    fflush(stderr);
#endif

    /* Symmetric Fincke-Pohst */
    for (i = 0; i < columns; i++) {
        fipo[i] = 0.0;
        dum1 = 0.0;
        for (j = 0; j < rows; j++) {
            tmp = 0.0;
            for (l = i; l < columns; l++)
                tmp += muinv[i][l]*bd[l][j]/c[l];
            fipo[i] += tmp*tmp;
            dum1 += fabs(tmp);
        }
        fipo[i] = SQRT(fipo[i]*Fd);
        dum1 =  fabs(dum1*Fq);
        if (dum1 < fipo[i])
            fipo[i] = dum1;
        fipo[i] *= (1.0+EPSILON);

        fipo_LB[columns][i] = -fipo[i];
        fipo_UB[columns][i] =  fipo[i];
#if VERBOSE > -1
        fprintf(stderr, "%0.3lf ", fipo[i]);
#endif
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
    for (i=columns-1; i>=0; i--) {
        if (fipo[i]<0.9) {
            printf("DEL\n");
            columns--;
        } else {
            break;
        }
    }
#endif
    /*|print_lattice();|*/
#if 0
    @<orthogonalize the basis@>;
    @<determine Fincke-Pohst bounds@>;
#endif

#if 0
    basis2LP(fipo_l,fipo_u);
#endif

    /* initialize first-nonzero arrays */
    for (l = 0; l < rows; l++) {
        for (i = 0; i < columns; i++) if (mpz_sgn(get_entry(i,l)) != 0) {
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
    hoelder_no = hoelder_success = 0;
    cs_success = nosolutions = loops = 0;
    N_success = 0;

    /* the loop of the exhaustive enumeration */
    do {
        /* increase loop counter */
        loops++;
        if ((stop_after_loops > 0) && (stop_after_loops <= loops))
            goto afterloop;

#if VERBOSE > -1
        if (loops % 100000000 ==0) {                 /*10000000*/
            printf("%ld loops, solutions: %ld ", loops, nosolutions);
#if FINCKEPOHST
            printf("fipo: %ld ", fipo_success);
#endif
            printf("\n");
            fflush(stdout);
        }
#endif

        /* compute new |cs| */
        olddum = dum[level];
        dum[level] = us[level] + y[level];
        cs[level] = cs[level+1] + dum[level]*dum[level]*c[level];

        if (isSideStep) {
            stepWidth = dum[level] - olddum;
        }

        if ((cs[level] < Fd) /*|&& (!prune0(fabs(dum[level]),N[level]))|*/)  {
#if FINCKEPOHST
            if (fabs(us[level]) > fipo[level]) {
                fipo_success++;
                goto side_step;
            }
#endif
            if (isSideStep) {
                compute_w2(w, bd, stepWidth, level, rows);
            } else {
                compute_w(w, bd, dum[level], level, rows);
            }

            if (level > 0) {
                /* not at a leave */
                i = prune_only_zeros(w, level, rows, Fq, first_nonzero_in_column, firstp,
                                     bd, y, us, columns);

                if (i < 0) {
                    goto step_back;
                } else if (i > 0) {
                    goto side_step;
                }

                if (prune(w[level], cs[level], rows, Fqeps)) {
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
                    isSideStep = 0;
                }
            } else {
                /* at $|level|=0$ */
                if (exacttest(w[0],rows,Fq) == 1) {
                    print_solution(w[level],rows,Fq,us,columns);
                    if ((stop_after_solutions > 0) && (stop_after_solutions <= nosolutions))
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
            isSideStep = 1;
        }
    } while (level<columns);
afterloop:

    /* final output */
    fprintf(stderr, "Prune_cs: %ld\n", cs_success);
    fprintf(stderr, "Prune_only_zeros: %ld of %ld\n", only_zeros_success, only_zeros_no);
    fprintf(stderr, "Prune_hoelder: %ld of %ld\n", hoelder_success, hoelder_no);
    fprintf(stderr, "Prune_N: %ld\n", N_success);
#if FINCKEPOHST
    printf("Fincke-Pohst: %ld\n", fipo_success);
#endif
    printf("Loops: %ld\n",loops);

    if ((stop_after_solutions <= nosolutions && stop_after_solutions > 0) ||
             (stop_after_loops <= loops && stop_after_loops > 0 )) {
        printf("Stopped after number of solutions: %ld\n", nosolutions);

        print_num_solutions(nosolutions);
        if ((stop_after_loops <= loops && stop_after_loops > 0)) {
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
    free (us);
    free (cs);
    free (y);
    free (delta);
    free (d);
    free (first_nonzero);
    free (first_nonzero_in_column);
    free (firstp);

    free (eta);
    free (v);
    for (l = 0; l <= columns; l++) free (w[l]);
    free (w);
    free (original_columns);

#if FINCKEPOHST
    free (fipo);
    for (l=0;l<columns;l++) free (muinv[l]);
    free(muinv);
#endif

    lllfree(mu,c,N,bd,columns);
    for (l=0;l<columns;l++) free (mu_trans[l]);
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
    cblas_daxpy(rows,alpha,bd[level],1,w[level],1);
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
    cblas_dcopy(rows,w[level+1],1,w[level],1);
    cblas_daxpy(rows,alpha,bd[level],1,w[level],1);
#else
    int l;

    l = rows - 1;
    while (l >= 0) {
        w[level][l] = w[level+1][l] + alpha*bd[level][l];
        l--;
    }
#endif
    return;
}

void gramschmidt(coeff_t **lattice, int columns, int rows, DOUBLE **mu,
                 DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq) {
    int i, l, j;
    DOUBLE dum;

    for (i = 0; i < columns; i++) {
        for (l = 0; l < rows; l++) bd[i][l] = (DOUBLE)mpz_get_d(get_entry(i,l));
        N[i] = 0.0;
        for (j=0;j<i;j++) {
            dum = 0.0;
            for(l=0;l<rows;l++) dum += (DOUBLE)mpz_get_d(get_entry(i,l)) * bd[j][l];
            mu[i][j] = dum / c[j];
            for (l = 0; l < rows; l++) bd[i][l] -= mu[i][j]*bd[j][l];
        }

        c[i] = scalarproductfp(bd[i],bd[i],rows);
        for(l = 0; l < rows; l++) N[i] += fabs(bd[i][l]);
        N[i] /= c[i];
        N[i] *= Fq;
#if VERBOSE > 0
        printf("%lf ",(double)c[i]);
#endif
/*|        N[i] *= (1.0 + EPSILON);|*/
    }
#if VERBOSE > 0
    printf("\n\n");
    fflush(stdout);
#endif
    return;
}

void givens(coeff_t **lattice, int columns, int rows, DOUBLE **mu,
            DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq) {
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
            mu[i][l] = (DOUBLE)mpz_get_d(get_entry(i,l));
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
        c[i] = mu[i][i]*mu[i][i];
        N[i] = 0.0;
        for (j = 0; j < rows; j++) {
            bd[i][j] *= mu[i][i];
            N[i] += fabs(bd[i][j]);
        }
        N[i] /= c[i];
        N[i] *= Fq;

        /* N[i] *= 1.0 + EPSILON; */

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
int exacttest(DOUBLE *v, int rows, DOUBLE Fq) {
    register int i;
    i = rows - 1;
    do {
        if (fabs(v[i]) > Fq*(1.0 + EPSILON)) {
            return 0;
        }
        i--;
    } while (i>=0);
    return 1;
}

int prune0(DOUBLE li, DOUBLE re) {
    if (li > re*(1+EPSILON)) {
        N_success++;
        return 1;
    } else {
        return 0;
    }
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

int print_solution(DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns) {
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

    if (!SILENT) {
        mpz_set_si(soltest_upfac,1);
        mpz_set_si(soltest_s,0);
        for (k=0;k<columns;k++) {
            if (ROUND(us[k])>0) {
                mpz_addmul_ui(soltest_s,get_entry(k,rows-1), ROUND(us[k]));
            } else {
            mpz_submul_ui(soltest_s,get_entry(k,rows-1),-ROUND(us[k]));
            }
        }

        i = 0;
        if (cut_after_coeff==-1) {
            end=no_original_columns;
        } else {
            end=cut_after_coeff;
        }
        for (j=0;j<end;j++) {
            if (original_columns[j]==0) {
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
                        mpz_addmul_ui(soltest_u,get_entry(k,i),ROUND(us[k]));
                    } else {
                        mpz_submul_ui(soltest_u,get_entry(k,i),-ROUND(us[k]));
                    }
                }
                mpz_sub(soltest_u,soltest_u,soltest_s);
                mpz_divexact(soltest_u,soltest_u,max_norm_initial);
                mpz_divexact(soltest_u,soltest_u,soltest_upfac);
                mpz_divexact_ui(soltest_u,soltest_u,denom);
                mpz_abs(soltest_u,soltest_u);
                if (!iszeroone && (mpz_cmp_si(soltest_u,0)<0 || mpz_cmp(soltest_u,upperbounds[i])>0) ) {  /* upperbounds not defined for $0/1$ problems */
                    fprintf(stderr," rounding error -> this is not a solutions!\n");
                    return 0;
                }
                i++;
            }
            mpz_out_str(NULL,10,soltest_u);
            fflush(stdout);
            mpz_out_str(fp,10,soltest_u);

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

void stopProgram(int sig) {
    if (sig != SIGALRM)
       return;

    printf("Stopped after SIGALRM, number of solutions: %ld\n", nosolutions);
    print_num_solutions(nosolutions);

    exit(11);
}

void shufflelattice() {
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
        for(i=lattice_columns-2; i>0; i--) {
            r = rand() % i;
            tmp = lattice[r];
            lattice[r] = lattice[i];
            lattice[i] = tmp;
        }
    }
    return;
}
