#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <gmp.h>
#include "dio2.h"
//#include "OpenBLAS/common.h"
//#include "OpenBLAS/cblas.h"

#define BLAS 1
#define USE_SSE 0
#define DEEPINSERT 1
#define DEEPINSERT_CONST 100
#define VERBOSE 1

#define GIVENS 1
#define LASTLINESFACTOR "1000000" /* "100000000" */
#define EPSILON 0.000001      /* 0.0001  */
#define LLLCONST_LOW  0.60 /* 0.75*/
#define LLLCONST_HIGH 0.90    /* 0.99 */
#define LLLCONST_HIGHER 0.999
#define ETACONST 0.55

#define SQRT sqrt
#define DOUBLE double
#define COEFF struct coe

/**
 * Definition of the lattice data structures
*/
struct coe {
    mpz_t c;
    int p;
};

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
mpz_t snd_q, snd_r, snd_s;

int system_rows, system_columns;
int lattice_rows, lattice_columns;
COEFF **lattice;
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
int iterate;
int no_iterates;
int bkz_beta, bkz_p;
int SILENT;
int nboundvars;

mpz_t soltest_u;
mpz_t soltest_s;
mpz_t soltest_upfac;

/**
 * Inline functions
 */
ROUND(r) ceil(r-0.5)
#define put_to(i,j,val) mpz_set(lattice[i][j+1].c,val)
#define smult_lattice(i,j,factor) mpz_mul(lattice[i][j+1].c,lattice[i][j+1].c,factor)
#define get_entry(i,j) lattice[i][j+1].c


@<lattice basis reduction algorithms@>;


long diophant(mpz_t **a_input, mpz_t *b_input, mpz_t *upperbounds_input,
    int no_columns, int no_rows,
    mpz_t factor_input, mpz_t norm_input, mpz_t scalelastlinefactor,
    int silent, int iterate, int iterate_no,
    int bkz_beta_input, int bkz_p_input,
    long stop_after_sol_input, long stop_after_loops_input,
    int free_RHS_input, int *org_col_input, int no_org_col_input,
    int cut_after, int nboundedvars, FILE* solfile) {

    int i,j;
    DOUBLE lD, lDnew;
    COEFF *swap_vec;

    /**
     * Initialize some globals
     */
    mpz_init_set(matrix_factor,factor_input);
    mpz_init_set(max_norm,norm_input);
    mpz_init(lastlines_factor);
    mpz_init(upfac);

    mpz_init(snd_q);
    mpz_init(snd_r);
    mpz_init(snd_s);

    mpz_init(soltest_u);
    mpz_init(soltest_s);
    mpz_init_set_ui(soltest_upfac,1);

    if (iterate) {
        no_iterates = iterate_no;
    } else {
        bkz_beta = bkz_beta_input;
        bkz_p = bkz_p_input;
    }
    SILENT = silent;
    stop_after_solutions = stop_after_sol_input;
    stop_after_loops = stop_after_loops_input;
    free_RHS = free_RHS_input;
    nom = 1;
    denom = 2;

    system_rows = no_rows;
    system_columns = no_columns;
    nboundvars = nboundedvars;

#if BLAS
    /*|goto_set_num_threads(1); |*/
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
    lattice = (COEFF**)calloc(lattice_columns,sizeof(COEFF*));
    for(j=0;j<lattice_columns;j++) {
        lattice[j] = (COEFF*)calloc(lattice_rows+1,sizeof(COEFF));
        for (i=0;i<=lattice_rows;i++) mpz_init(lattice[j][i].c);
    }

    /**
     * read the system
     */
    for (j=0;j<system_rows;j++) {
        for (i=0;i<system_columns;i++) {
            mpz_mul(lattice[i][j+1].c,a_input[j][i],matrix_factor);
        }
        mpz_mul(lattice[system_columns][j+1].c,b_input[j],matrix_factor);
    }

    /**
     * handle upper bounds
     */
    mpz_init_set_si(upperbounds_max,1);
    iszeroone = 1;
    if (upperbounds_input==NULL) {
        fprintf(stderr, "No upper bounds: 0/1 variables are assumed \n"); fflush(stderr);
    } else {
        upperbounds = (mpz_t*)calloc(system_columns,sizeof(mpz_t));
        for (i=0;i<system_columns;i++) mpz_init_set_si(upperbounds[i],1);
        for (i=0;i<nboundvars/*|system_columns|*/;i++) {
            mpz_set(upperbounds[i],upperbounds_input[i]);
            if (mpz_sgn(upperbounds[i])!=0) {
                mpz_lcm(upperbounds_max,upperbounds_max,upperbounds[i]);
            }
        }
        if (mpz_cmp_si(upperbounds_max,1)>0) iszeroone = 0;
        fprintf(stderr,"upper bounds found. Max="); fflush(stderr);
        mpz_out_str(stderr,10,upperbounds_max);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    /**
     * handle preselected columns
     */
    if (org_col_input!=NULL)
        no_original_columns = no_org_col_input;
    else
        no_original_columns = system_columns;

    original_columns = (int*)calloc(no_original_columns,sizeof(int));

    if (org_col_input!=NULL)
        for (i=0;i<no_original_columns;i++) original_columns[i] = org_col_input[i];
    else {
        for (i=0;i<no_original_columns;i++) original_columns[i] = 1;
        fprintf(stderr, "No preselected columns \n"); fflush(stderr);
    }

    /**
     * append the other parts of lattice
     */
    for (j=system_rows;j<lattice_rows;j++) {
        mpz_mul_si(lattice[j-system_rows][j+1].c,max_norm,denom);
        mpz_mul_si(lattice[lattice_columns-2][j+1].c,max_norm,nom);
    }
    mpz_set(lattice[system_columns+free_RHS][lattice_rows].c,max_norm);

    if (free_RHS) {
        mpz_set_si(lattice[system_columns][lattice_rows-1].c,1);
        mpz_set_si(lattice[system_columns+1][lattice_rows-1].c,0);
    }
    mpz_set(lattice[system_columns+free_RHS][lattice_rows].c,max_norm);
    for (i=0;i<lattice_columns-1;i++) coeffinit(lattice[i],lattice_rows);

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
    for (i=lattice_columns-2;i>0;i--) lattice[i] = lattice[i-1];
    lattice[0] = swap_vec;

#if 0
    printf("After permute\n");
    print_lattice();
#endif
    shufflelattice();
    /**
     * first reduction
     */
    mpz_set_ui(lastlines_factor,1);
    fprintf(stderr, "\n"); fflush(stderr);
    lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_LOW);

#if 0
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

#if 0
    printf("After cutting\n");
    print_lattice();
#endif

#if 1
    shufflelattice();
    /**
     * second reduction
     */
    mpz_set_ui(lastlines_factor,1);
    lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_HIGH);
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
    /*|mpz_set_str(lastlines_factor,LASTLINESFACTOR,10);|*/
    mpz_set(lastlines_factor, scalelastlinefactor);
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
    if (iterate) {
        iteratedlll(lattice,lattice_columns-1,lattice_rows,no_iterates,LLLCONST_HIGH);
    } else {
        shufflelattice();
        lDnew = bkz(lattice,lattice_columns,lattice_rows,LLLCONST_HIGHER,40,bkz_p);

        i = 0;
        do {
            lD = lDnew;
            shufflelattice();
            lDnew = bkz(lattice,lattice_columns,lattice_rows,LLLCONST_HIGH,bkz_beta,bkz_p);
            printf("%0.3lf %0.3lf %0.3lf\n",lD, lDnew, lD-lDnew);
            i++;
        }
        while (i<1 && fabs(lDnew-lD)>0.01);
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

#if 0
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

 void stopProgram() {
     printf(stderr, "Stopped after SIGALRM, number of solutions: %ld\n", nosolutions);
     print_num_solutions(nosolutions);

     exit(11);
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

void coeffinit(COEFF *v, int z) {
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

int lllfp (COEFF **b, DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs,
            int start, int s, int z, DOUBLE delta) {

    int i, j, k;
    DOUBLE ss;
    int counter;

    int Fc, Fr;
    DOUBLE mus, cc;
    mpz_t musvl;
    mpz_t hv;
    DOUBLE *swapd;

    int ii, iii;
    COEFF *swapvl;
    COEFF *bb;

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

    counter = 0;

    /* The main loop */
    while (k<s) {
#if VERBOSE > 3
        if ((counter%500)==0) @+ {
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
        if (k == 1)
            c[0] = N[0] * N[0];
        c[k] = N[k] * N[k];
        for (j = 0; j < k; j++) {
            ss = scalarproductfp(bs[k],bs[j],z);
            if (fabs(ss) < N[k] * N[j] / TWOTAUHALF) {
                ss = (DOUBLE)scalarproductlfp(b[k],b[j]);
            }
            for (i = 0; i < j; i++)
                ss -= mu[j][i] * mu[k][i] * c[i];
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
                mpz_set_d(musvl,mus);
                if (fabs(mus)>TWOTAUHALF) {
                    Fc = 1;
                }

                Fr = 1;
                /* set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$ */
                switch (mpz_get_si(musvl)) {
                case 1:
                    /* $\lceil\mu_{k,j}\rfloor = 1$ */
                    i = b[j][0].p;
                    while (i != 0) {
                            bb = &(b[k][i]);
                            mpz_sub(bb->c, bb->c, b[j][i].c);
                            iii = bb->p;
                            if ((b[k][i-1].p != i) && (mpz_sgn(bb->c) != 0))
                                for (ii= i - 1; (ii >=0) && (b[k][ii].p == iii); ii--)
                                    b[k][ii].p = i;
                            else if (mpz_sgn(bb->c) == 0) {
                                for (ii = i - 1;  (ii >= 0) && (b[k][ii].p == i); ii--)
                                    b[k][ii].p = iii;
                            }
                            i = b[j][i].p;
                    }
                    for(i=0;i<j;i++) mu[k][i] -= mu[j][i];
                    break;

                case -1:
                    /* $\lceil\mu_{k,j}\rfloor = -1$ */
                    i = b[j][0].p; @;
                    while (i != 0) {
                            bb = &(b[k][i]);
                            mpz_add(bb->c, bb->c, b[j][i].c);
                            iii = bb->p;
                            if ((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0))
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==iii);ii--) b[k][ii].p = i;
                            else if (mpz_sgn(bb->c)==0) {
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==i);ii--) b[k][ii].p = iii;
                            }
                            i = b[j][i].p;
                    }
                    for(i=0;i<j;i++) mu[k][i] += mu[j][i];
                    break;

                default:
                    /* $\lceil\mu_{k,j}\rfloor \neq \pm 1$ */
                    i=b[j][0].p; @;
                    while (i!=0) {
                            bb=&(b[k][i]);
                            mpz_submul(bb->c,b[j][i].c,musvl);
                            iii = bb->p;
                            if ((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0))
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==iii);ii--) b[k][ii].p = i;
                            else if (mpz_sgn(bb->c)==0) {
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==i);ii--) b[k][ii].p = iii;
                            }
                            i = b[j][i].p;
                    }
#if 0
                    daxpy(j,-mus,mu[k],1,mu[j],1);
#endif
                    for(i=0;i<j;i++) mu[k][i] -= mu[j][i]*mus;

                }


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
        @<fourth step: deepinsert columns@>;
#else
        @<fourth step: swap columns@>;
#endif
    }
    mpz_clear(hv);
    mpz_clear(musvl);
    return(1);

}
