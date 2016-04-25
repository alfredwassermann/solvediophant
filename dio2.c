#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <gmp.h>
#include "diophant.h"
#include "OpenBLAS/common.h"
#include "OpenBLAS/cblas.h"

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

/**
 * Inline functions
 */
ROUND(r) @[ceil(r-0.5)@]
#define put_to(i,j,val) @[mpz_set(lattice[i][j+1].c,val)@]
#define smult_lattice(i,j,factor) @[mpz_mul(lattice[i][j+1].c,lattice[i][j+1].c,factor)@]
#define get_entry(i,j) @[lattice[i][j+1].c@]


@<basic subroutines@>;
@<lattice basis reduction algorithms@>;


long diophant(mpz_t **a_input, mpz_t *b_input, mpz_t *upperbounds_input,
    int no_columns, int no_rows,
    mpz_t factor_input, mpz_t norm_input, mpz_t scalelastlinefactor,
    int silent, int iterate, int iterate_no,
    int bkz_beta_input, int bkz_p_input,
    long stop_after_sol_input, long stop_after_loops_input,
    int free_RHS_input, int *org_col_input, int no_org_col_input,
    int cut_after, int nboundedvars, FILE* solfile)@;
{
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
     * set the lattice dimension@>;
     */
    lattice_rows = system_rows + system_columns + 1;
    lattice_columns = system_columns + 2;

    if (free_RHS) {
        lattice_rows++;
        lattice_columns++;
    } else {
        fprintf(stderr,"The RHS is fixed !\n");@+ fflush(stderr);
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
#ifndef NO_OUTPUT
        printf("No upper bounds: 0/1 variables are assumed \n"); @+ fflush(stdout);
#endif
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
#ifndef NO_OUTPUT
        fprintf(stderr,"upper bounds found. Max=");
        mpz_out_str(stderr,10,upperbounds_max);
        fprintf(stderr,"\n");
        @+ fflush(stderr);
#endif
    }

    /**
     * handle preselected columns
     */
    if (org_col_input!=NULL) no_original_columns = no_org_col_input;
    else no_original_columns = system_columns;

    original_columns = (int*)calloc(no_original_columns,sizeof(int));

    if (org_col_input!=NULL)
        for (i=0;i<no_original_columns;i++) original_columns[i] = org_col_input[i];
    else {
        for (i=0;i<no_original_columns;i++) original_columns[i] = 1;
        printf("No preselected columns \n"); fflush(stdout);
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
    mpz_init_set_si(max_up,1); @;
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

    @<undo scaling of last rows@>;
    for (i = 0; i < lattice_columns; i++)
        mpz_divexact(lattice[i][lattice_rows].c,lattice[i][lattice_rows].c,lastlines_factor);
    if (free_RHS)
        for (i = 0; i < lattice_columns; i++)
            mpz_divexact(lattice[i][lattice_rows-1].c,
                lattice[i][lattice_rows-1].c,lastlines_factor);
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
     * close solution file@>;
     */
    if (SILENT) fprintf(fp, "%ld solutions\n", nosolutions);
    fflush(fp);

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
