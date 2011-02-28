\datethis
@*diophant -- Solving Diophantine Linear Systems.
The library {\tt diophant} enables the solution of 
Diophantine Linear Systems based on LLL reduction.
The user has to provide an integer matrix $A$, a right hand side
vector $b$, optionally a vector $u$ of upper bounds on the solution vector $x$.
{\tt diophant} computes all integer solutions $x$, such that
$$A\cdot x = b, \hbox{ where } 0\leq x\leq u.$$
An example is the following $3\times 13$-system, to be solved with a
0/1-vector:
$$
\left(\matrix{
5&5&5&5&5&1&1&1&0&0&0&0&0\cr
0&4&4&8&0&0&0&0&4&4&4&0&0\cr
2&1&4&3&2&1&1&1&2&5&2&3&1\cr}\right) \cdot x = 
\left(\matrix{12\cr 12\cr 12\cr}\right)
$$
One 0/1-solution is 
$$x= (1,0,0,1,0,1,0,1,1,0,0,1,0)^\top.$$

@ Initial definitions. 
The LLL and BKZ algorithms can be controlled by the following
declarations.
@d DEEPINSERT 1
@d DEEPINSERT_CONST 2000
@d VERBOSE 1

@d GIVENS 1
@d LASTLINESFACTOR "1" /* "100000000" */
@d EPSILON 0.000001      /* 0.0001  */
@d LLLCONST_LOW  0.75 /* 0.75*/
@d LLLCONST_HIGH 0.99    /* 0.99 */
@d LLLCONST_HIGHER 0.999 

@d SQRT sqrt
@d DOUBLE double
@d COEFF struct coe

@f mpz_t long
@f DOUBLE double
@f COEFF int

@c
@<include header files@>;
@<definition of the lattice data structures@>;
@<global variables@>;
@<inline functions@>;
@<basic subroutines@>;
@<lattice basis reduction algorithms@>;

@ The main routine {\tt diophant()}: 
@c
long diophant(mpz_t **a_input, mpz_t *b_input, mpz_t *upperbounds_input,
    int no_columns, int no_rows,
    mpz_t factor_input, mpz_t norm_input, 
    int silent, int iterate, int iterate_no,
    int bkz_beta_input, int bkz_p_input, 
    long stop_after_sol_input, long stop_after_loops_input, 
    int free_RHS_input, int *org_col_input, int no_org_col_input,
    int cut_after, int nboundedvars, FILE* solfile)@;
{
    int i,j;
    DOUBLE lD, lDnew;
    COEFF *swap_vec;
    
    @<initialize some globals@>;
    @<set the lattice dimension@>;
    @<allocate memory@>;
    @<read the system@>;
    @<handle upper bounds@>;
    @<handle preselected columns@>;
    @<append the other parts of lattice@>;
    @<open solution file@>;
#if 0
    printf("Before scaling\n");
    print_lattice();
#endif
    @<scale lattice@>;
#if 0
    printf("After scaling\n");
    print_lattice();
#endif
    
#if 0    
    print_NTL_lattice();   /* Version for the NTL output */
    return 0;
#endif

    @<permute lattice columns@>;
#if 0
    printf("After permute\n");
    print_lattice();
#endif
    shufflelattice();
    @<first reduction@>;
#if 0
    printf("After first reduction\n");
    print_lattice();
#endif
    @<cut the lattice@>;
#if 0
    printf("After cutting\n");
    print_lattice();
#endif
#if 1
    shufflelattice();
    @<second reduction@>;
#endif
#if 0
    printf("After second reduction\n");
    print_lattice();
#endif
    
    @<scale last rows@>;
    @<third reduction@>;
    @<undo scaling of last rows@>;
    
#if 0
    printf("Before enumeration\n");
    print_NTL_lattice();   /* Version for the NTL output */
    /*|print_lattice();|*/
#endif

    @<explicit enumeration@>;
    @<close solution file@>;
    @<free multiprecision memory@>;
    return nosolutions;
}

@ The header file {\tt diophant.h}. It is needed if one wants to include 
|diophant()| in his program.
@(diophant.h@>=
#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h>
#undef BLAS
extern long diophant(mpz_t **a_input, mpz_t *b_input, mpz_t *upperbounds_input,
    int no_columns, int no_rows,
    mpz_t factor_input, mpz_t norm_input, 
    int silent, int iterate, int iterate_no,
    int bkz_beta_input, int bkz_p_input, 
    long stop_after_sol_input, long stop_after_loops_input, 
    int free_RHS_input, int *org_col_input, int no_org_col_input,
    int cut_after, int nboundedvars, FILE* solfile);
#endif

@ The following header files are included. There is a conflict
between freelip and the BLAS-library. So if we use multiprecision
integers we cannot use BLAS.
@<include header files@>=
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <gmp.h>
#undef BLAS
#if defined(BLAS)
#include "blas/blas1.h"
#else
#if defined(BLAS2)
#include "blas_intel/cblas.h"
#endif
#endif

@ The data structure of the lattice:
@<definition of the lattice data structures@>=
struct coe {
    mpz_t c;
    int p;
};

@ The global variables.
We start with the 
variables which control the scaling of the lattice.
|max_norm| resembles the $\lambda$ in design problems.
@<global variables@>=
    mpz_t matrix_factor;      
    mpz_t max_norm;
    mpz_t max_norm_initial;
    mpz_t max_up;
    mpz_t dummy;
    long nom, denom;
    mpz_t lastlines_factor;

@ The variables which define the lattice and its dimensions.
@<global variables@>=
    int system_rows, system_columns;
    int lattice_rows, lattice_columns;
    COEFF **lattice;
    int free_RHS;
    int iszeroone;
    mpz_t *upperbounds;
    mpz_t upperbounds_max;
    mpz_t upfac;

@ Other variables.
@<global variables@>=
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
    
@ The access to the lattice is done via inline functions.
The rounding function is defined as
@d ROUND(r) @[ceil(r-0.5)@]

@<inline functions@>=
#define put_to(i,j,val) @[mpz_set(lattice[i][j+1].c,val)@]
#define smult_lattice(i,j,factor) @[mpz_mul(lattice[i][j+1].c,lattice[i][j+1].c,factor)@]
#define get_entry(i,j) @[lattice[i][j+1].c@]

@ Read the parameters of {\tt diophant} and initialize some
global variables.
@<initialize some globals@>=
    mpz_init_set(matrix_factor,factor_input);
    mpz_init_set(max_norm,norm_input);
    mpz_init(lastlines_factor);
    mpz_init(upfac);

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
    
@ @<free multiprecision memory@>=
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

@ The lattice has two more columns than the original system.
The number of rows is the number of columns plus number of rows
plus one of the original system. 
If the right hand side is free (e.g. for design problems)
there is an additional column and an additional row.
@<set the lattice dimension@>=
    lattice_rows = system_rows + system_columns + 1;
    lattice_columns = system_columns + 2;
    
    if (free_RHS) {
        lattice_rows++;
        lattice_columns++;
    } else {
#ifndef NO_OUTPUT
        fprintf(stderr,"The RHS is fixed !\n");@+ fflush(stderr);
#endif
    }
    cut_after_coeff = cut_after;
    
@ Allocate memory for the lattice array and fill it with zero.
@<allocate memory@>=
    lattice = (COEFF**)calloc(lattice_columns,sizeof(COEFF*));
    for(j=0;j<lattice_columns;j++) {
        lattice[j] = (COEFF*)calloc(lattice_rows+1,sizeof(COEFF));
        for (i=0;i<=lattice_rows;i++) mpz_init(lattice[j][i].c); 
    }
@ The input matrix and the right hand side
  vector are copied row by row into the lattice.
@<read the system@>=
    for (j=0;j<system_rows;j++) {
        for (i=0;i<system_columns;i++) {
            mpz_mul(lattice[i][j+1].c,a_input[j][i],matrix_factor);
        }
        mpz_mul(lattice[system_columns][j+1].c,b_input[j],matrix_factor);
    }
@ The upper bounds are copied. If the input vector points to |NULL|, a 
0/1 system is assumed. The variable 
|upperbounds_max| contains the least common multiple
of the upper bounds, let's call it $u_{\max}$. 
Starting from $u_{\max}=1$ 
it is computed with $0\le i<n$:
$$
u_{\max} \leftarrow u_{\max}\cdot {u_i\over \gcd(u_{\max},u_i)}.
$$
It is needed for an
appropiate scaling of the lattice, see section |@<scale lattice@>|.
Further, this determines if we have a 0/1 system, i.~e. if
$|upperbounds_max|=1$.
@<handle upper bounds@>=
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
@ Handle the original columns. This is used for preselected columns.
If the array points to |NULL|, we fill the array completely with 1's.
This is of course only useful in case of 0/1 problems.
The non-0/1 case has still to be handled.
@<handle preselected columns@>=
    if (org_col_input!=NULL) no_original_columns = no_org_col_input;
    else no_original_columns = system_columns;

    original_columns = (int*)calloc(no_original_columns,sizeof(int));

    if (org_col_input!=NULL) 
        for (i=0;i<no_original_columns;i++) original_columns[i] = org_col_input[i];
    else {
        for (i=0;i<no_original_columns;i++) original_columns[i] = 1;
#ifndef NO_OUTPUT
        printf("No preselected columns \n"); @+ fflush(stdout); 
#endif
    }

@ The matrix of unity, the last columns and rows are appended.
The diagonal (unity) matrix is set to
$|denom|*|max_norm|$ which is $2\cdot\lambda$ in design problems
and $2$ in problems with fixed right hand side.
The second last column (the column corresponding to the right hand side)
is set to
$|nom|*|max_norm|$ which is $\lambda$ in design problems
and $1$ in problems with fixed right hand side.
@<append the other parts of lattice@>=
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
    
@ The lower parts of the lattice are multiplied by factors stemming
from the upper bounds. This is only necessary if we have non-0/1 bounds.
Each row in the lower part of the lattice corresponds to a variable.
We multiply the diogonal entry in row $j$ (of the lower part,
therefore in the whole system it is row  $j+ |system_rows|$) by 
${u_{\max}\over u_j}$.
The right hand side column entry in row $j$ is multiplied by
$u_{\max}$.
@<scale lattice@>=
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
    
@ The lattice columns are cyclically shifted: 
the second last column becomes the first column,
the first column becomes the second, $\ldots$
@<permute lattice columns@>=
    swap_vec = lattice[lattice_columns-2];
    for (i=lattice_columns-2;i>0;i--) lattice[i] = lattice[i-1];
    lattice[0] = swap_vec;
    
@ The first reduction: Pure LLL to compute
    the kernel. 
@<first reduction@>=
    mpz_set_ui(lastlines_factor,1);
#ifndef NO_OUTPUT
    printf("\n"); @+ fflush(stdout);
#endif
    lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_LOW); 

@ The second reduction: Pure LLL with higher quality.
This is done in order to find solutions by chance.
@<second reduction@>=
    mpz_set_ui(lastlines_factor,1);
    lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_HIGH); 
#ifndef NO_OUTPUT
    printf("Second reduction successful\n"); @+ fflush(stdout);
#endif

@ Delete the unnecessary parts of the reduced lattice, multiply the last
    row(s) by an appropiate factor, such that
    there is only one entry in these rows after the next reduction.
    This is done after first lattice basis reduction.
    If it was not successful in computing a complete basis of the kernel
    we leave |diophant()|. 
    The function |cutlattice()| is defined in section
    |@<cut lattice@>|.
@<cut the lattice@>=    
    if (cutlattice()) {
#ifndef NO_OUTPUT
        printf("First reduction successful\n"); @+ fflush(stdout);
#endif
    } else {
#ifndef NO_OUTPUT
        printf("First reduction not successful\n"); @+ fflush(stdout);
#endif
        return 0;
    }
    
    for (j=0;j<lattice_columns-1;j++) solutiontest(j);
@ The last row or in case of a free right hand side the two last rows
are multiplied by a large factor |LASTLINESFACTOR|. Then
it easier to search with the exhaustive enumeration for
nonhomogeneous solutions because the columns corresponding to the 
right hand side (and the additional column in case of a free right hand side)
appear only once in the basis vectors.
{bf Attention:} It only works for fixed right hand side

@<scale last rows@>=
    mpz_set_str(lastlines_factor,LASTLINESFACTOR,10);
    for (i=0;i<lattice_columns;i++) 
        mpz_mul(lattice[i][lattice_rows].c,lattice[i][lattice_rows].c,lastlines_factor);
    if (free_RHS) 
        for (i=0;i<lattice_columns;i++) 
            mpz_mul(lattice[i][lattice_rows-1].c,lattice[i][lattice_rows-1].c,lastlines_factor);

@ Undo the scaling of section |@<scale last rows@>|
after the second reduction.
@<undo scaling of last rows@>=
    for (i=0;i<lattice_columns;i++) 
        mpz_divexact(lattice[i][lattice_rows].c,lattice[i][lattice_rows].c,lastlines_factor);
    if (free_RHS) 
        for (i=0;i<lattice_columns;i++) 
            mpz_divexact(lattice[i][lattice_rows-1].c,
                lattice[i][lattice_rows-1].c,lastlines_factor);

@ The third reduction is done with blockwise Korkine Zolotareff reduction.
The last column of |lattice| is a zero vector.
@<third reduction@>=
#ifndef NO_OUTPUT
    printf("\n"); @+ fflush(stdout);
#endif
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
#ifndef NO_OUTPUT
    printf("Third reduction successful\n"); @+ fflush(stdout);
#endif

@ The exhaustive enumeration of all solutions.
@<explicit enumeration@>=
#ifndef NO_OUTPUT
    printf("\n"); @+ fflush(stdout);
#endif
    nosolutions = explicit_enumeration(lattice,lattice_columns-1,lattice_rows);

@*Subroutines for {\tt diophant}.
@<basic subroutines@>=
@<debug print@>;
@<print the lattice@>;
@<shuffle lattice columns@>;
@<print NTL lattice@>;
@<gcd@>;
@<Initialize the lattice@>;
@<cut lattice@>;
@<solutiontest@>;

@ @<debug print@>=
void debug_print(char *m, int l) {
#ifndef NO_OUTPUT
    if (VERBOSE>=l) {
        printf("debug>> %s\n",m); @+ fflush(stdout);
    }
#endif
    return;        
}

@ Print the lattice for debugging purposes.
The lattice is printed columnwise (i.e.~in transposed form).
@<print the lattice@>=
#if 1
void print_lattice() {    
    int i,j;
    for (i=0;i<lattice_columns;i++) {
        for (j=0;j<lattice_rows;j++) {
            mpz_out_str(NULL,10,get_entry(i,j));
            printf(" ");
        }
        printf("\n");
    }
    printf("\n"); @+ fflush(stdout);
    return;
}
#else
void print_lattice() {    
    int i,j;
    for (j=0;j<lattice_rows;j++) {
        for (i=0;i<lattice_columns;i++) {
            mpz_out_str(NULL,10,get_entry(i,j));
            printf(" ");
        }
        printf("\n");
    }
    printf("\n"); @+ fflush(stdout);
    return;
}
#endif

@ Print the lattice in NTL format.
The lattice is printed columnwise (i.e.~in transposed form).
@<print NTL lattice@>=
void print_NTL_lattice() {  
    int i,j;
#if 1
    printf("%d\n",system_rows);
    printf("\n[");
    for (i=0;i<lattice_columns-1;i++) {   /* Don't print the last vector (only zeroes). */
        printf("[");
        for (j=0;j<lattice_rows;j++) {
            mpz_out_str(NULL,10,get_entry(i,j));
            printf(" ");
        }
        printf("]");
        printf("\n");
    }
    printf("]\n"); @+ fflush(stdout);
#if 0    
    printf("\n");
    mpz_out_str(NULL,10,upperbounds_max);
    printf("\n\n[");
    for (i=0;i<lattice_columns-2;i++) {
        mpz_out_str(NULL,10,upperbounds[i]);
        printf(" ");
    }
    printf("]\n"); @+ fflush(stdout);
#endif    
#endif    
    return;
}

@ Compute GCD in case of single precision arithmetics.
This piece of code was written by Brendan McKay.
@<gcd@>=
long gcd(long n1,long n2) {
        long a,b,c;

        if (n1 > n2) {
            a = n1; @+ b = n2;
        } else {
            a = n2; @+ b = n1;
        }

        while ((c = a % b) > 0) {
            a = b; @+ b = c;
        }
        return b;
}
@ Construct the information of the sparse structure of the lattice.
@<Initialize the lattice@>=
void coeffinit(COEFF *v, int z)
{
    short r = 0;
    short i;
    for (i=z;i>=0;i--) {
        v[i].p = r;
        if (mpz_sgn(v[i].c)!=0) r=i;
    }
    return;
}

@ Delete the unnecessary columns and rows of the lattice after
the first reduction, i.e.~only the basis of the kernel remains.
Returns $0$ if it is not possible to achieve a solution.
Finally, the unnecessary rows are deleted and the sparse structure of the
lattice is updated.
@<cut lattice@>=
int cutlattice() {
    int j,i,flag;

    @<delete unnecessary columns@>;
    @<test for right hand side columns@>;
    
    for (j=0;j<lattice_columns;j++)  {     /* Now the rows are deleted. */
       if (nboundvars==0) {
            for (i=system_rows;i<lattice_rows;i++) 
                put_to(j,i-system_rows,get_entry(j,i));
        } else {
            for (i=system_rows;i<system_rows+nboundvars;i++) 
                put_to(j,i-system_rows,get_entry(j,i));
            for (i=system_rows+system_columns;i<lattice_rows;i++) 
                put_to(j,i-system_rows-system_columns+nboundvars,get_entry(j,i));
        }
    }
    lattice_rows -= system_rows;
    lattice_rows -= (system_columns-nboundvars);

     for (j=0;j<lattice_columns;j++) coeffinit(lattice[j],lattice_rows);
    
    return 1;
}

@ Shuffle the coluns of the lattice.
@<shuffle lattice columns@>=
void shufflelattice() {
    COEFF *tmp;
    int i, j, r;
    
    srand((unsigned)(time(0))); 
    for (j=0;j<100;j++) {
        for(i=lattice_columns-2; i>0; i--) {
            r = rand() % i;
            tmp = lattice[r];
            lattice[r] = lattice[i];
            lattice[i] = tmp;
        }
    }
    return;
}

@ Delete the columns which do not belong to the kernel.
@<delete unnecessary columns@>=
    j=0;
    do {
        if (lattice[j][0].p>system_rows)
            j++;
        else {
            for (i=j+1;i<lattice_columns;i++) lattice[i-1] = lattice[i];
            lattice_columns--;
        }
    } while (j<lattice_columns-1);

@ Test if the remaining columns contain a nonzero entry
in the last row. Otherwise a nonhomogenous solution is not possible.

@<test for right hand side columns@>=
    flag = 0;
    for (i=0;i<lattice_columns;i++) if (mpz_sgn(get_entry(i,lattice_rows-1))!=0) {
        flag = 1;
        break;
    }
    if (flag==0) {
#ifndef NO_OUTPUT
        printf("Nonhomogenous solution not possible.\n"); @+ fflush(stdout);             
#endif
        exit(1);        
        return 0;
    }            
@ Test of column |position| for a solution during the reduction phase.
@<solutiontest@>=
int solutiontest(int position) {
    int i,j;
    int low, up;    
    int end;

    @<test the last two rows@>;
    @<test, if column is a solution@>;

        mpz_set_si(upfac,1);
        mpz_divexact(soltest_s,get_entry(position,lattice_rows-1),lastlines_factor);
    @<write a solution with blanks@>;
    @<test if one solution is enough@>;
    return 1;
}
@ @<global variables@>=
    mpz_t soltest_u;
    mpz_t soltest_s;
    mpz_t soltest_upfac;
@ @<initialize some globals@>=
    mpz_init(soltest_u);
    mpz_init(soltest_s);
    mpz_init_set_ui(soltest_upfac,1);

@ Test the last two rows. 
@<test the last two rows@>=
    if (mpz_cmpabs(get_entry(position,lattice_rows-1),max_norm)!=0) return 0;
    if (mpz_sgn(get_entry(position,lattice_rows-1-free_RHS))==0) return 0;
@ A final test, if the column is really a solution.
We have to distinguish whether this subroutine is called during the first or
the secoond reduction. 
    If $|lattice_columns|= |system_columns|+2+|free_RHS|$
    is true, then no columns have been deleted so far, i.~e.~ we are
    in the first reduction phase. Accordingly a start index
    |low| is determined.
@<test, if column is a solution@>=
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
@ Print a solution. The variables are separated with blanks.
Only in the case that just one solution is needed, the solution is also written
into the solution file.
@<write a solution with blanks@>=
    i = low;
    if (cut_after_coeff==-1) {
        end=no_original_columns;
#if 0        
        if (nboundvars!=0) { /* conflicts with |original_columns| */
            end=nboundvars;
        }
#endif        
    } else {
        end=cut_after_coeff;
    }
        
    for (j=0;j<end;j++) {
        if (original_columns[j]==0) {
            mpz_set_si(soltest_u,0);
        } else {
            if (!iszeroone) {
                if (mpz_cmp_si(upperbounds[i-low],0)!=0) {
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
    printf("\n"); @+ fflush(stdout);
    if (stop_after_solutions==1) fprintf(fp,"\n"); 

@ The case $|stop_after_solutions|=1$: We can stop already.
@<test if one solution is enough@>=
    if (stop_after_solutions==1) {
#ifndef NO_OUTPUT
        printf("Stopped in phase 1 after finding a random solution\n");
#endif
        exit(0);
    }
@*Elementary lattice basis reduction. We have 
standard lattice basis reduction (with deep insertions) and blockwise
Korkine Zolotareff reduction. The main ingredient is the subroutine 
|lllfp| which does lattice basis reduction.
@<lattice basis reduction algorithms@>=
@<lll-subroutines@>;
@<the underlying lattice basis reduction@>;
@<compute $\log D$@>;
@<orthogonal defect@>;
@<standard-lll@>;
@<iterated-lll@>;
@<bkz-reduction@>;
@<overall exhaustive enumeration@>;

@ We begin with the core |lllfp|.
The multiprecision version needs additionally 
the matrix |bs|. This matrix contains the floating point
approximation of the lattice |b|.
@<the underlying lattice basis reduction@>=
#define TWOTAUHALF 67108864.0 /* $2^{\tau/2}$*/
int lllfp (COEFF **b, DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, 
            int start, int s, int z, DOUBLE delta)@;
{
    @<variables for |lllfp|@>; @#
    mpz_init(musvl);
    mpz_init(hv);

    if ((z <= 1) || (s <= 1)) {    /* Test for trivial cases. */
#ifndef NO_OUTPUT
        printf("Wrong dimensions in lllfp\n"); @+ fflush(stdout);
#endif
        return(0);
    }
    
    k = (start > 1) ? start : 1; @#
    
    @<first step: compute norms@>; @#
    counter = 0;
    while (k<s) {   /* The main loop. */
#if VERBOSE > 3
        if ((counter%500)==0) @+ {
#ifndef NO_OUTPUT
            @+printf("LLL: %d k:%d\n",counter,k); @+ fflush(stdout);
#endif
        } 
        counter++;
#endif
        @<second step: orthogonalization of $b_k$ @>;
        @<third step: size reduction of $b_k$@>;
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
@ Step 1:
We begin with stage $k=1$. Schnorr's original algorithm runs
from $k=2$ to $k=s$. Here we have $k=1,\ldots,s-1$.

Initially the $b_i$'s are rounded and then their norms,
$N_i = ||b'_i||$, are computed.

Then we do steps 2 -- 4 while  $k\leq s$.
@<first step: compute norms@>=
    if (k<1) k = 1;
    for (i=k-1;i<s;++i) {
        ss = 0.0;
         for (j=0;j<z;++j) {
            bs[i][j] = (DOUBLE)mpz_get_d(b[i][j+1].c);
            ss += bs[i][j]*bs[i][j];
        }
        N[i] = SQRT(ss);
    }
@ @<variables for |lllfp|@>=
    int      i,j,k;
    DOUBLE   ss;
    int      counter;
@ Step 2:
computation of $\mu_{k,1},\ldots,\mu_{k,k-1}$ and $c_k = ||\hat{b}_k||^2.$
This is done with Gram Schmidt Orthogonalization.
@<second step: orthogonalization of $b_k$@>=
        if (k==1) c[0] = N[0]*N[0];
        c[k] = N[k]*N[k];
        for (j=0;j<k;j++) {
            ss = scalarproductfp(bs[k],bs[j],z);
            if (fabs(ss)< N[k]*N[j]/TWOTAUHALF) {
                ss = (DOUBLE)scalarproductlfp(b[k],b[j]);
            } 
            for (i=0;i<j;i++) ss -= mu[j][i]*mu[k][i]*c[i];
            mu[k][j] = ss/c[j];
            c[k] -= ss*mu[k][j];
        } 
@ Step 3: size reduction. This one of the two main parts of
LLL reduction. It is some kind of rounded Gram Schmidt orthogonalization
of the integer basis of the lattice.

The flag |Fr| is set true if the column has been changed, i.~e.~if there really
has been a size reduction.
The flag |Fc| is set true in |@<round the Gram Schmidt coefficient@>| if the
Gram-Schmidt coefficient is too large.

At the end of this section we test for 
rounding errors and linear dependencies.
If $F_c$ is
true, then the stage is decreased by one. Otherwise
the new value of $b'_k$ is computed and we proceed to step 4.
@<third step: size reduction of $b_k$@>=
        Fc = Fr = 0;         
        for(j=k-1;j>=0;j--) {
            if (fabs(mu[k][j])>0.5) {
                @<round the Gram Schmidt coefficient@>;
                Fr = 1;  
                @<set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$@>;
                mu[k][j] -= mus;
                  solutiontest(k);
            }   
        } 
        @<recompute $N_k$@>;
        if (Fc==1) {
            k = (k-1 > 1) ? k-1 : 1;    /* $k=\max(k-1,1)$ */
        } else {
            @<test for linear dependencies@>;
        }
                
@ @<variables for |lllfp|@>=
    int      Fc, Fr;
    DOUBLE   mus, cc;
    mpz_t musvl;
    mpz_t hv;
    DOUBLE   *swapd; 

@ The Gram-Schmidt coefficient is rounded. If it is too large (|Fc| is true),
we have to do Schnorr's correction step:
We will step back due to rounding errors.
@<round the Gram Schmidt coefficient@>=
        mus = ROUND(mu[k][j]);
                mpz_set_d(musvl,mus);
                if (fabs(mus)>TWOTAUHALF) {
                    Fc = 1;      
#if 0                    
                    printf("correct possible rounding errors\n");@+fflush(stdout);
#endif                    
                }
                    /* We have to correct possible rounding errors.*/

@ @<variables for |lllfp|@>=
    int ii, iii;
    COEFF    *swapvl;
    COEFF        *bb;
@ Replace $b_k$ by $b_k \leftarrow b_k - \mu b_j$ and replace
   $\mu_{k,i}$ by $\mu_{k,i} \leftarrow \mu_{j,i}\mu$,
    where $\mu = \lceil\mu_{k,j}\rfloor$ as computed in section
    |@<round the Gram Schmidt coefficient@>|.
   This computation uses the sparse structure of the vectors.
   
   In order to save multiplications, we treat the three cases 
   $\lceil\mu_{k,j}\rfloor = \pm 1$ and $\lceil\mu_{k,j}\rfloor \neq \pm 1$
    separately.
@<set $b_k = b_k - \lceil\mu_k,j\rfloor b_j$@>=
                switch (mpz_get_si(musvl)) {
                case 1:
                    @<$\lceil\mu_{k,j}\rfloor = 1$@>;
                    break;
                    
                case -1:
                    @<$\lceil\mu_{k,j}\rfloor = -1$@>;
                    break;
                    
                default:
                    @<$\lceil\mu_{k,j}\rfloor \neq \pm 1$@>;
                }
@ The original loop was
    |for(i=1;i<=z;i++) b[k][i].c -=  b[j][i].c;|
@<$\lceil\mu_{k,j}\rfloor = 1$@>=
                    i=b[j][0].p; 
                    while (i!=0) { 
                            bb=&(b[k][i]); 
                            mpz_sub(bb->c,bb->c,b[j][i].c);
                            iii = bb->p;
                            if ((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0))
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==iii);ii--) b[k][ii].p = i;
                            else if (mpz_sgn(bb->c)==0) {
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==i);ii--) b[k][ii].p = iii;
                            }
                            i = b[j][i].p;
                    }
                    for(i=0;i<j;i++) mu[k][i] -= mu[j][i];
@ The original loop was
    |for(i=1;i<=z;i++) b[k][i].c +=  b[j][i].c;|
@<$\lceil\mu_{k,j}\rfloor = -1$@>=
                    i=b[j][0].p; @;
                    while (i!=0) {
                            bb=&(b[k][i]);
                            mpz_add(bb->c,bb->c,b[j][i].c);
                            iii = bb->p;
                            if ((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0)) 
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==iii);ii--) b[k][ii].p = i;
                            else if (mpz_sgn(bb->c)==0) {
                                for (ii=i-1;(ii>=0)&&(b[k][ii].p==i);ii--) b[k][ii].p = iii;
                            }
                            i = b[j][i].p;
                    }
                    for(i=0;i<j;i++) mu[k][i] += mu[j][i];
@ The original loop was 
|for(i=1;i<=z;i++) b[k][i].c -=  b[j][i].c*musvl;|
@<$\lceil\mu_{k,j}\rfloor \neq \pm 1$@>=
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
@ @<recompute $N_k$@>=
    {
            N[k] = 0.0;
            for(i=0;i<z;i++) {
                bs[k][i] = (DOUBLE)mpz_get_d(b[k][i+1].c);
                N[k] += bs[k][i]*bs[k][i];
            }   
            N[k] = SQRT(N[k]);
    }
@ Before going to step 4 we test if $b_k$ is linear dependent.
This is the case if $||N_k||<{1\over 2}$.
If we found a linear dependent vector $b_k$, 
we shift $b_k$ to the last column of the
matrix and restart |lllfp| with $s:= s-1$.
@<test for linear dependencies@>=
            if (N[k]<-EPSILON) {
#ifndef NO_OUTPUT
                fprintf(stderr,"Nk negativ! contact the author.\n"); 
                @+ fflush(stderr); 
                printf("Nk negativ! contact the author.\n"); @+ fflush(stdout); @+
#endif
                exit(1);
            }
            if (N[k]<0.5) {
               swapvl = b[k];
               ss = N[k];
               swapd = bs[k];
               for(i=k+1;i<s;i++) {
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

@ @<variables for |lllfp|@>=
    int Fi;

@ Step 4: Swapping of columns with deep insertions.
Will be used if the compiler flag |DEEPINSERT| is set.
The vector $b_k$ will be inserted left as possible.
Swapping of the matrix columns is done entirely with 
pointers arithmetics. 
If no insertion is possible
we step forward to
$k\leftarrow k+1$, otherwise we go back to $k\leftarrow k-1$.
@<fourth step: deepinsert columns@>=
            cc = N[k]*N[k];
            j = 0;
            Fi = 0;
            while (j<k) {
#if 0
                if ((j>DEEPINSERT_CONST && j<k-DEEPINSERT_CONST) || delta*c[j]<=cc) { 
#endif                
                if (delta*c[j]<=cc) { 
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
            if (Fi==1) k = (j-1 > 1) ? j-1 : 1;       /* $k = \max(j-1,1)$ */
            else {
                k++;
            }

@ Standard exchange of columns
(original Step 4):
Either swap the vectors $b_{k-1}$ and  $b_k$ or increment $k$.
Swapping of the matrix columns is done entirely with swapping of
pointers. 
If the $b_k$ and $b_{k-1}$ are exchanged, we step back to
$k\leftarrow k-1$, otherwise we go forward to $k\leftarrow k+1$.
@<fourth step: swap columns@>=
            if (delta*c[k-1] > c[k]+mu[k][k-1]*mu[k][k-1]*c[k-1]) {
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
            } else k++;

@*Subroutines needed for lattice basis reduction.
These are scalarproducts, norms, allocation and orthogonalization.
@<lll-subroutines@>=
@<scalarproduct with integers@>;
@<scalarproduct with doubles@>;
@<memory allocation@>;
@<free the memory@>;

@ The integer scalarproduct is able to use the sparse structure of the
   vectors. But it strongly depends on the hardware and the compiler.
   This is the new variant which uses the sparse structure of the input 
vectors.

It depends on the compiler and the machine whether this loop
is fast enough to gain speed against the plain scalarproduct evaluation.
$t_1$ runs through the vector $v$, $t_2$ runs through the vector $w$.
@<scalarproduct with integers@>=
DOUBLE scalarproductlfp (COEFF *v, COEFF *w )
{
    DOUBLE    erg;
    long    t1,t2;
    COEFF    *vv,*ww;

    erg = 0.0;
    t1 = v[0].p;
    t2 = w[0].p;
    if ((t1==0)||(t2==0)) return 0;

    do {
        if (t2>t1) {
            t1 = v[t2-1].p;
            if (t2!=t1) { 
                if (t1==0) break;
                t2 = w[t2].p; 
                if (t2==0) break;
            }
            else goto gleich;
        }
        else if (t2<t1) {
            t2 = w[t1-1].p;
            if (t2!=t1) { 
                if (t2==0) break;
                t1 = v[t1].p; 
                if (t1==0) break;
            }
            else goto gleich;
        }
        else {
gleich:    vv = &(v[t1]);
            ww = &(w[t2]);
            erg += (DOUBLE)mpz_get_d(vv->c) * (DOUBLE)mpz_get_d(ww->c);
            t1 = vv->p;
            if (t1==0) break;
            t2 = ww->p;
            if (t2==0) break;
        }
    }
    while (1);

    return (erg);
}
@ @<scalarproduct with doubles@>=
DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w , int n)
{
    DOUBLE r;
#if defined(BLAS)
    r = ddot(n,v,1,w,1);
#else
#if defined(BLAS2)
    r = cblas_ddot(n,v,1,w,1);
#else
    int i;

    r = 0.0;
    for (i=n-1;i>=0;i--) r += v[i]*w[i];
#endif    
#endif    
    return r;
}

@ Allocation of memory. For the upper triangular matrix $\mu$ a
rectangular matrix is allocated. So it fits for orthogonalization with
Gram-Schmidt, Householder and Givens rotation.
@<memory allocation@>=
int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N,  DOUBLE ***bs, int s, int z) @;
{
    int i,m;

    if ((z < 1) || (s < 1)) return 0;

    (*c)=(DOUBLE*)calloc(s,sizeof(DOUBLE));
    (*N)=(DOUBLE*)calloc(s,sizeof(DOUBLE));
    (*mu)=(DOUBLE**)calloc(s,sizeof(DOUBLE*));
    for(i=0;i<s;i++) (*mu)[i]=(DOUBLE*)calloc(z,sizeof(DOUBLE));
    
    m = (z>s) ? z : s;
    (*bs)=(DOUBLE**)calloc(m,sizeof(DOUBLE*));
        
    for(i=0;i<m;i++) 
        (*bs)[i]=(DOUBLE*)calloc(m,sizeof(DOUBLE));

    return 1;
}
@ @<free the memory@>=
int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s) @;
{
    int i;

    for (i=0;i<s;++i) free(bs[i]);
    free(bs);
    for (i=0;i<s;++i) free(mu[i]);
    free(mu);
    free(N);
    free(c);

    return 1;
}

@ Compute $\log D$ for a lattice $L$, see the LLL paper.
@<compute $\log D$@>=
double logD(COEFF **lattice, DOUBLE *c, int s, int z) {
    double d = 0.0;
    int i;
    for (i=0;i<s;i++) {
        d += log(c[i])*(s-i);
    }
    d *= 0.5;
    return d;
}

@ The logarithm of the orthogonality defect of a lattice is computed:
$d = {\prod_{i=0}^{s-1} \Vert b_i\Vert\over \det(L)}$.
For $\det(L)$ we use $\prod_{i=0}^{s-1} \Vert\hat b_i\Vert$.
@<orthogonal defect@>=
double orthogonal_defect(COEFF **lattice, DOUBLE *c, int s, int z) {
    double defect = 0.0;
    
#if 0
    int i;
    for (i=0;i<s;i++) defect += log((double)normfp(lattice[i])) 
        - log((double)c[i]);
#endif    
    defect /= 2.0;
    return defect;
}
    
@ Standard LLL reduction
@<standard-lll@>=
void lll (COEFF **b, int s, int z, DOUBLE quality)
{
    DOUBLE **mu;
    DOUBLE *c;
    DOUBLE *N;
    DOUBLE **bs;
    int r;

    lllalloc(&mu,&c,&N,&bs,s,z);
    r = lllfp(b,mu,c,N,bs,1,s,z,quality);
    lllfree(mu,c,N,bs,s);

    return;
}

@ Iterated LLL reduction. Standard LLL-reduction is applied to
the lattice. Then the columns are sorted in descending (or ascending) order and
the process is iterated. The euclidean length of the lattice 
vectors is stored in |N|.
@<iterated-lll@>=
DOUBLE iteratedlll(COEFF **b, int s, int z, int no_iterates, DOUBLE quality)
{
    DOUBLE **mu;
    DOUBLE *c;
    DOUBLE *N;
    DOUBLE **bs;
    int r,l,i,j, runs;
    COEFF *swapvl;
    DOUBLE lD;

    lllalloc(&mu,&c,&N,&bs,s,z);
    r = lllfp(b,mu,c,N,bs,1,s,z,quality);

    lD = logD(b,c,s,z);
#ifndef NO_OUTPUT
    printf("   log(D)= %f\n", lD); 
    @+ fflush(stdout);
#endif

    for (runs=1;runs<no_iterates;runs++) {


       for (j=s-1;j>0;j--) {
           for (l=j-1;l>=0;l--) { 
            @|
                /*|if (N[l]<N[j]) {|*/    /* $<$ sorts 'in descending order.' */
                if (N[l]>N[j]) {    /* $>$ sorts 'in ascending order.' */
                     swapvl = b[l];
                    for (i=l+1;i<=j;i++) b[i-1] = b[i];
                    b[j] = swapvl;
                }
            }
        }
    
        r = lllfp(b,mu,c,N,bs,1,s,z,quality);
        lD = logD(b,c,s,z);
#ifndef NO_OUTPUT
        printf("%d: log(D)= %f\n",runs, lD); 
        @+ fflush(stdout);
#endif
    }
    
    lllfree(mu,c,N,bs,s);

    return lD;
}

@*Blockwise Korkine Zolotareff reduction.
@<bkz-reduction@>=
@<compute gamma function@>;
@<exhaustive enumeration of a block@>;
@<bkz-algorithm@>;

@ The control algorithm of blockwise Korkine Zolotarev Reduction.
There are 2 parameters $\beta$ and $\delta$ with ${1\over 4}<\delta<1$ (same
as in LLL) and $2<\beta<m$ (the blocksize).
The parameter |s| is the number of columns of the lattice,
where the last column is the zero vector.
Therefore the last vector of the mathematical lattice is
$|last|\leftarrow |s|-2$.
|bkz| uses the subroutines |lllfp| and |enumerate|.
$p$ controls the pruning in the pruned Gauss enumeration.
$|start_block|$ is cyclically shifted through $0,1,\ldots,s-2$.
   |zaehler| counts the number of positions $j$ that satisfy
   $\delta||\hat{b}_j||^2 \leq \lambda_1(\pi_j(L(b_j,\ldots,b_k)))$.
   |zaehler| is reset to $0$ if the inequality does not hold for $j$.

 Step 1:
 $\mu$ and $c$ have to be allocated with $s+1$ columns.
@<bkz-algorithm@>=
DOUBLE bkz(COEFF **b, int s, int z, DOUBLE delta, int beta, int p)
{
     @<bkz variables@>;

     last = s-2;    /* |last| points to the last nonzero vector of the lattice.*/
     if (last<1) {
#ifndef NO_OUTPUT
        printf("BKZ: the number of basis vectors is too small.\n");
        printf("Probably the number of rows is less or equal");
        printf(" to number of columns in the original system\n");
        printf("Maybe you have to increase c0 (the first parameter)!\n");
#endif
        mpz_clear(hv);
        return 0.0;
    }

    u = (long*)calloc(s,sizeof(long));
    for (i=0;i<s;i++) u[i] = 0;

    lllalloc(&mu,&c,&N,&bs,s,z);
    lllfp(b,mu,c,N,bs,1,s,z,delta);

    start_block = zaehler = -1;
    while (zaehler<last) {
    
    start_block++;
    if (start_block==last) start_block = 0;
    end_block = (start_block+beta-1 < last) ? start_block+beta-1 : last; 
        /*  $|end_block| := \min(|start_block|+\beta-1,|last|)$  */
#if 0
printf("start_block=%d, end_block=%d\n",start_block,end_block);
#endif
    new_cj = enumerate(mu,c,u,s,start_block,end_block,p);
        /* The exhaustive enumeration. */
    h = (end_block+1 < last) ? end_block+1 : last; 
        /*  $h := \min(|end_block|+1,|last|)$  */

    if (delta*c[start_block] > new_cj) {
        @<successful enumeration@>;
        zaehler = -1;
    }
    else {
        if (h>0) {
            lllfp(b,mu,c,N,bs,h-2,h+1,z,delta);   /* For some unkown reason we have to
                                                 use $h-2$ as |start|. */
        }
        zaehler++;
    }
  } /* end of |while| */

    lD = logD(b,c,s-1,z);
#ifndef NO_OUTPUT
    printf("bkz: log(D)= %f\n", lD); 
    @+ fflush(stdout);
#endif
     lllfree(mu,c,N,bs,s);
     free(u);
     mpz_clear(hv);
     return lD;
}
@ @<bkz variables@>=
    DOUBLE **mu, *c, *N;
    DOUBLE **bs; 
    static mpz_t hv;
    int zaehler;
    int h,i,last;
    int start_block, end_block;
    long *u;
    DOUBLE new_cj;
    DOUBLE lD;
    
    mpz_init(hv);

@ Now we found a new vector whose corresponding Gram Schmidt vector
   is shorter than $c_j$ by at least a factor $\delta$.
    Write the new linear combination at position $j-1$ after shifting all 
   vectors behind $j-1$ one position to the right.

   If $N_{h-1} < 0.5$  the vectors in $b$ were linear dependent.
   Since at most one of the vectors $b_j,\ldots,b_{h}$ was linear dependent
   and is now zero at position $h$, swap the vectors back.
   Otherwise the vectors weren't linear dependent. This is the case if a
   multiple of $b_{s-1}$ was found in |enumerate|.
@<successful enumeration@>=
#if defined(ORIGINAL_SCHNORR_EUCHNER)
    @<old integration of the new vector, part 1@>;
#else
    @<build new basis@>; 
#endif
@#
        lllfp(b,mu,c,N,bs,start_block-1,h+1,z,delta);

        if (N[h]<-EPSILON) {
#ifndef NO_OUTPUT
            fprintf(stderr,"NN negativ\n"); @+ fflush(stderr); 
            printf("NN negativ\n"); @+ fflush(stdout); @+
#endif
            exit(1);
        }
#if defined(ORIGINAL_SCHNORR_EUCHNER) 
    @<old integration of the new vector, part 2@>;
#endif

@ @<bkz variables@>=
    int g, ui, q, j;
    COEFF *swapvl;

@ The new vector 
$b_a = \sum_{i=a}^e u_ib_i$ is computed, where
$a= |start_block|$ and $e= |end_block|$.
The other old vectors $b_{a+1},\ldots,b_e$ are adapted,
in order to keep the whole block linear independent.
The algorithm is described in the doctoral thesis of
H.~Ritter.
@<build new basis@>=
    for (j=1;j<=z;j++) mpz_set_si(b[last+1][j].c,0);
    for (i=start_block;i<=end_block;i++) {    
        if (u[i]!=0) for (j=1;j<=z;j++) {
            if (u[i]>0) {
                mpz_addmul_ui(b[last+1][j].c,b[i][j].c,u[i]);
            } else {
                mpz_submul_ui(b[last+1][j].c,b[i][j].c,-u[i]);
            }
        }
    }
    g = end_block;
    while (u[g]==0) g--;
    
    i = g-1;
    while (labs(u[g])>1) {
        while (u[i]==0) i--;
        q = (int)ROUND((1.0*u[g])/u[i]);
        ui = u[i];
        u[i] = u[g] - q*u[i];
        u[g] = ui;
                
        for (j=1;j<=z;j++) {
            mpz_set(hv,b[g][j].c);
            mpz_mul_si(b[g][j].c,b[g][j].c,(long)q);
            mpz_add(b[g][j].c,b[g][j].c,b[i][j].c);
            mpz_set(b[i][j].c,hv);
        }
        coeffinit(b[g],z);
        coeffinit(b[i],z);
    }
    
    swapvl = b[g];
    for (i=g;i>start_block;i--) b[i] = b[i-1];
    b[start_block] = b[last+1];
    coeffinit(b[start_block],z);
    
    b[last+1] = swapvl;
    for (j=1;j<=z;j++) mpz_set_si(b[last+1][j].c,0);
    coeffinit(b[last+1],z);
    
@ This is the original way to  include the new vector:
  just write it to the beginning of the block and reduce
  the new system. Then one vector is reduced to zero.
  This vector is deleted.
@<old integration of the new vector, part 1@>=
        for (l=1;l<=z;l++) mpz_set_si(b[last+1][l].c,0);
        for (i=start_block;i<=end_block;i++) 
            for (l=1;l<=z;l++) {
                mpz_addmul_si(b[last+1][l].c,b[i][l].c,ui);
            }                
        coeffinit(b[last+1],z);
        solutiontest(last+1);
        swapvl = b[last+1];
        for (i=last;i>=start_block;i--) b[i+1] = b[i];
        b[start_block] = swapvl;

@ @<old integration of the new vector, part 2@>=
        if (N[h]<0.5) {     /* After reduction this vector should be zero */
            swapvl = b[h];
            for (i=h+1;i<=last+1;i++) b[i-1] = b[i];
            b[last+1] = swapvl;
        } else {
#ifndef NO_OUTPUT
            printf("Not linear dependent; %f\n",(double)(N[h-1])); @+ fflush(stdout); 
#endif
            exit(1);
        }        

@ Pruned Gauss-Enumeration.
Enumerate in depth first search. Pruned Version without goto's 
as it is described in H.~H.~H\"orner's diploma thesis.
Blocks with blocksize lower than |SCHNITT| are completely enumerated,
otherwise the enumeration is pruned.
$2^{-p}$ is the probability that lattice points are lost. 
@<exhaustive enumeration of a block@>=
DOUBLE enumerate(DOUBLE **mu, DOUBLE *c, long *u, int s, @|
int start_block, int end_block, int p)
{
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
#ifndef NO_OUTPUT
          fprintf(stderr,"Hier ist was faul! start_block=%d %f\n",start_block,
        (double)c[start_block]); @+ fflush(stderr); 
          printf("Hier ist was faul! start_block=%d %f\n",start_block,
        (double)c[start_block]); @+ fflush(stdout);
#endif        
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
     for (i=0;i<s;i++) {
        sigma[i] =(DOUBLE*)calloc(s,sizeof(DOUBLE));
        r[i] = i-1;
     }
         
@#
     len = (end_block+1-start_block);
     for (i=start_block;i<=end_block+1;i++) {
         cs[i] = y[i] = 0.0;
         u[i] = us[i] = v[i] = delta[i] = 0;
         d[i] = 1;
     }
     us[start_block] = u[start_block] = 1;
     cd = c[start_block];

     t = tmax = start_block;      /* Now we start from $t=|start_block|$ instead of 
      $t=|end_block|$. */

     @<precompute $\eta$@>;
     while (t<=end_block) {
        @<the block search loop@>;
     } 
     free (us);
     free (cs);
     free (y);
     free (delta);
     free (d);
     free (eta);
     free (v);
     for (i=s-1;i>=0;i--) {
        free(sigma[i]);
     }
     free(sigma);
     free(r);
     return (cd);
}

@ The search loop.
@<the block search loop@>=
{
        dum = us[t] + y[t];
        cs[t] = cs[t+1] + dum*dum*c[t];

#if 0     
        if (cs[t]<cd-eta[t]) {
        /*|alpha = (t>t2) ? 0.7*cd : cd;|*/
#else        
        if (len <= SCHNITT) {
            alpha = 1.0;
        } else {
            alpha = sqrt(1.20*(end_block+1-t)/len);
            if (alpha>=1.0) alpha = 1.0;
        }
        alpha *= cd;
        /*|
        alpha = (2*t>len) ? 0.8*cd : cd;
        |*/
        if (cs[t] < alpha-EPSILON) {
#endif            
            if (t>start_block) {
                t--;
                if (r[t+1]>r[t]) r[t] = r[t+1];
                
                delta[t] = 0;
                for (i=r[t+1];i>t;i--) sigma[i][t] = sigma[i+1][t] + us[i]*mu[i][t];
#if 0                
                dum = 0.0;
                for (i=t+1;i<=tmax;i++) dum += us[i]*mu[i][t];
                if (fabs(dum-sigma[t+1][t])>0.001) {
                    printf("1diff: %0.6lf %0.6lf %0.8lf\n", dum, sigma[t+1][t], dum-sigma[t+1][t]);
                }
#endif                
                y[t] = sigma[t+1][t]; /*|dum;|*/
                us[t] = v[t] = (long)(ROUND(-y[t]));
                d[t] = (v[t] > -y[t]) ? -1 : 1;
            } else {
printf("success %0.4lf %0.4lf %0.8lf\n",cd, cs[start_block], cd-cs[start_block]);
                cd = cs[start_block];
                for (i=start_block;i<=end_block;i++) u[i] = us[i];
                goto nextstep;                
            }
        } else {
            t++;
            r[t] = t;
nextstep: @;
            if (tmax<t) tmax = t;
            if (t<tmax) delta[t] = -delta[t];
            if (delta[t]*d[t]>=0) delta[t] += d[t];
            us[t] = v[t] + delta[t];
        }
}          

@ Precomputation of $\eta_t$. 
If we know $c$, an upper bound of the square of the norm 
of a shortest vector, and
we arrived in our enumeration at level $t$ and a vector whose square norm
equals $\tilde{c}_t$, then the probability that we this branch of the
enumeration tree yields a shorter vector than the bound $c$
is
${{\rm vol}\,S(\sqrt{c-\tilde{c}_t},z)\over \det L }
$.
We compute $\eta_t$ such that 
$$
{{\rm vol}S(\sqrt{\eta},z)\over \det L } = 2^{-p}.
$$
Then, in our enumeration we stop if
$$
    \tilde{c}_t\ge c -\eta_t.
$$
The volume of the $z$-dimensional ball with radius $\sqrt{\eta}$
is equal to
$$
\hbox{vol}\,S(\sqrt{\eta},z)=
{(\pi\cdot\eta)^{z/2}\over \Gamma(z/2+1)}
$$ 
The dimension $z$ is $t-1$, one less than the enumeration level, and
we have $t=i-|start_block|$. 
H.H.~H\"orner computes $\eta$ by an approximation 
with Stirlings's formular. Here it is computed with the %'
Euler-MacLaurin-formular, see
Abramowitz, Stegun.

Update 19.1.2010:
now we try H\"orners approximation:
$$
\eta_i = {i-j\over 2\pi e}\biggl(\pi(i-j)p^2\sum_{h=j}^i c_h\biggr)^{1/(i-j)}
$$

@<precompute $\eta$@>=
    eta[start_block] = 0.0;
     if (end_block-start_block <= SCHNITT) {
         for (i=start_block+1;i<=end_block;i++) eta[i] = 0.0;
     }
     else {
         dum = log(c[start_block]);
/* This is my version to cokkmpute the parameters of the Gaussian volume
   heuristics. The Hoerner version seems to be better.
*/
#if 0         
         for (i=start_block+1;i<=end_block;i++) {
            t_up = i-start_block;
             eta[i] = exp( (log_gamma(t_up/2.0+1)-p*log(2.0))*2.0/t_up + dum/t_up )/pi;
             if (i<end_block) dum += log(c[i]);
#if 0
            if (i==start_block+1) printf("eta: \n");
            printf("%0.2lf ",eta[i]);
            if (i==end_block) printf("\n");
#endif
         }
#else
/*
    Hoerners version of the Gaussian volumne heuristics.
*/
        dum = log(c[start_block]);
        for (i=start_block+1;i<=end_block;i++) {
            t_up = i-start_block;
            eta[i] = 0.5*t_up*exp( (log(pi*t_up)-2.0*p*log(2.0)+dum)/t_up )/(pi*e);
            if (i<end_block) dum += log(c[i]);
#if 0
            if (i==start_block+1) printf("eta: ");
            printf("%0.2lf ",eta[i]);
            if (i==end_block) { 
                printf("\n");
                fflush(stdout);
            }
#endif
        }
#endif        
     }
@ @<compute gamma function@>=
DOUBLE laurin(DOUBLE x) {
    static DOUBLE K1 =0.9181938533204672741780329736405620;
    static DOUBLE K2 =0.0833333333333333333333333333333333333;
    static DOUBLE K3 =0.0027777777777777777777777777777777777;
    static DOUBLE K4 =0.000793650793650793650793650793650793650;
    static DOUBLE K5 =0.0005952380952380952380952380952380952380;
    static DOUBLE K6 =0.000841750841750841750841750841750841750;
    static DOUBLE K7 =0.001917526917526917526917526917526917526;
    static DOUBLE K8 =0.00641025641025641025641025641025641025;
    static DOUBLE K9 =0.0295506529510021209716796875;
    static DOUBLE K10 =0.17968122661113739013671875;
    static DOUBLE K11 =1.39243221282958984375;
    DOUBLE y;
    
    y = 1.0/(x*x);
    y = (x-0.5)*log(x) - x + K1 + 
   (1.0/x)*(K2-y*(K3-y*(K4-y*(K5-y*(K6-y*(K7-y*(K8-y*(K9-y*(K10-y*K11)))))))));
    return y;
}

DOUBLE log_gamma(DOUBLE x) {
    DOUBLE y;
    int i,n;
    static int MM = 13;
    
    if (x<=0.0) return -1.0;
    
    if (x>100000000.0) { 
        y = x*(log(x)-1.0);
    } else {
        if (x>=MM) {
            y = laurin(x);
        } else {
            n = MM - (int)(floor(x));
            y = x - floor(x) + MM;
            y = laurin(y);
            for (i=0;i<n;i++) y -= log(y+i);
        }    
    }
    return y;
}

@*Exhaustive enumeration. 
The algorithm of H.~Ritter.
@d FINCKEPOHST 1
@d FINCKEPOHSTEXTREME 0

@<overall exhaustive enumeration@>=
    @<globals for enumeration@>;
    @<some vector computations@>;
    @<orthogonalization@>;
    @<matrix inversion for Fincke-Pohst@>;
    @<pruning subroutines and output@>;
    @<output polyhedra@>;
    @<output LP@>;
    
DOUBLE explicit_enumeration (COEFF **lattice, int columns, int rows)
{
    @<local variables for |explicit_enumeration()|@>;
    
    long nlow[1000];
    for (i=0;i<1000;i++) nlow[i] = 0;


    @<test the size of the basis@>;
    @<allocate the memory for enumeration@>;
    @<initialize arrays@>;
    @<count nonzero entries in the last rows(s)@>;
#if 1
    @<sort lattice columns@>;
#endif    
    @<set the simple pruning bounds@>;
   @<orthogonalize the basis@>;
#if 0    
    print_lattice();
    basis2poly();
#endif    
    
#if defined(FINCKEPOHST)
    @<determine Fincke-Pohst bounds@>;
#if FINCKEPOHSTEXTREME
    @<sort by Fincke-Pohst bounds@>;
    @<orthogonalize the basis@>;
    @<determine Fincke-Pohst bounds@>;
#endif
#endif

#if 0
    basis2LP(fipo_l,fipo_u);
#endif

    @<initialize first-nonzero arrays@>;
    @<more initialization@>;
    @<the loop of the exhaustive enumeration@>;

    @<final output@>;
    @<free allocated memory for enumeration@>;
    return 1;
}
    
@ @<globals for enumeration@>=
#if 0
static FILE *fp;
#endif 

@ @<local variables for |explicit_enumeration()|@>=
    int level,level_max;
    int i,j,l;
    long loops;
    
    DOUBLE *y, *cs, *us;

    long *delta, *d, *eta;
#if 0        
        mpz_t *v;
#else        
        long *v;
#endif        
    int *first_nonzero, *first_nonzero_in_column, *firstp; 

    DOUBLE *N, **mu, *c, **w, **bd;

    DOUBLE Fd, Fq; 
    DOUBLE dum;
    COEFF *swap_vec;

    DOUBLE **sigma;
    int *r;
    
#if defined(FINCKEPOHST)
    DOUBLE *fipo, *fipo_u, *fipo_l;
#endif

@ It is tested, if the remaining columns build a basis of the kernel.
@<test the size of the basis@>=
#ifndef NO_OUTPUT
    printf("Dimension of solution space (k): %d compared to s-z+2: %d\n",
        columns, system_columns-system_rows+1+free_RHS);
    fflush(stdout);
#endif    
    
    if (columns < system_columns-system_rows+1+free_RHS) {
#ifndef NO_OUTPUT
        fprintf(stderr,"LLL didn't succeed in computing a basis of the kernel.\n");
        fprintf(stderr,"Please increase c0 (the first parameter)!\n");
        printf("LLL didn't succeed in computing a basis of the kernel.\n");
        printf("Please increase c0 (the first parameter)!\n");
#endif
        return 0;
    }
    

@  The memory is allocated. If multiprecision arithmetics is used
|bd| has already been allocated in |lllalloc|.
@<allocate the memory for enumeration@>=
    lllalloc(&mu,&c,&N,&bd,columns,rows); 

    us=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    cs=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    y=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    delta=(long*)calloc(columns+1,sizeof(long));
    d=(long*)calloc(columns+1,sizeof(long));
    first_nonzero=(int*)calloc(rows,sizeof(int));
    first_nonzero_in_column=(int*)calloc(columns+rows+1,sizeof(int));
    if (first_nonzero_in_column == NULL) return(0);
    firstp=(int*)calloc(columns+1,sizeof(int));

    eta=(long*)calloc(columns+1,sizeof(long));
    v=(long*)calloc(columns+1,sizeof(long));
    w=(DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
    for (i=0;i<=columns;i++) w[i]=(DOUBLE*)calloc(rows,sizeof(DOUBLE));

    sigma =(DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
    r =(int*)calloc(columns+1+1,sizeof(int));
    for (i=0;i<columns+1;i++) {
       sigma[i] =(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
       r[i] = i-1;
    }

#if defined(FINCKEPOHST)
    fipo=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    fipo_u=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    fipo_l=(DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
    muinv=(DOUBLE**)calloc(columns,sizeof(DOUBLE*));
    for(i=0;i<columns;++i) muinv[i]=(DOUBLE*)calloc(rows,sizeof(DOUBLE));
#endif
#if 0
    upb =(mpz_t*)calloc(columns+1,sizeof(mpz_t));
    lowb =(mpz_t*)calloc(columns+1,sizeof(mpz_t));
    for (i=0;i<=columns;i++) {
        mpz_init(upb[i];
        mpz_init(lowb[i];
    }
#endif    

@ @<globals for enumeration@>=
#if defined(FINCKEPOHST)
    DOUBLE **muinv;    
#endif    
    /*|mpz_t *upb,*lowb;|*/
    long fipo_success;
    
@    The starting positions of the arrays are set.
@<initialize arrays@>=
    for (i=0;i<=columns;i++) {
        cs[i] = y[i] = us[i] = 0.0;
        delta[i] = 0;
#if 0                
                mpz_set_si(v[i],0);
#else
                v[i] = 0;
#endif                
        eta[i] = d[i] = 1;
        for (l=0;l<rows;l++) w[i][l] = 0.0;
    }

@ Do simple sorting along the last row such that the 
    columns with a nonzero entry in the last row are 
    at the end.    
@<sort lattice columns@>=
   for (j=columns-1;j>0;j--) {
    for (l=j-1;l>=0;l--) {
        if (mpz_cmpabs(get_entry(l,rows-1),get_entry(j,rows-1))>0) {
                swap_vec = lattice[l];
                for (i=l+1;i<=j;i++) lattice[i-1] = lattice[i];
                lattice[j] = swap_vec;
            }
        }
    }

@ Do simple sorting along Fincke-Pohst bounds.
@<sort by Fincke-Pohst bounds@>=
   for (j=columns-3;j>0;j--) {
       for (l=j-1;l>=0;l--) {
#if 0        
        if (fipo[l]<fipo[j]) {      /* $<$ means: sort in descending order. */
#else
        if (fipo_u[l]-fipo_l[l]<fipo_u[j]-fipo_l[j]) {      
#endif            
            swap_vec = lattice[l];
            for (i=l+1;i<=j;i++) lattice[i-1] = lattice[i];
            lattice[j] = swap_vec;
                }
        }
    }

@ @<count nonzero entries in the last rows(s)@>=
    if (free_RHS) {
        i=0;
        for (j=columns-1;j>=0;j--) if (mpz_sgn(get_entry(j,rows-2))!=0) i++;
#ifndef NO_OUTPUT
        printf("Number of nonzero entries in the second last row: %d\n",i);
        fflush(stdout);
#endif
    }

    i=0;
    for (j=columns-1;j>=0;j--) if (mpz_sgn(get_entry(j,rows-1))!=0) i++;
#ifndef NO_OUTPUT
    printf("Number of nonzero entries in the last row: %d\n",i);
    fflush(stdout);
#endif

@    Set the simple bounds for pruning.
    |Fq| is the maximum norm of a solution,
    |Fd| is the square of the euclidean norm of a solution,
    |EPSILON| is a security buffer to avoid pruning
    due to rounding errors.
@<set the simple pruning bounds@>=
    Fq = (DOUBLE)mpz_get_d(max_norm); 
    Fd = (rows*Fq*Fq)*(1.0+EPSILON); 
#ifndef NO_OUTPUT
#if VERBOSE > 0
    printf("Fq: %f\n",(double)Fq);
    printf("Fd: %f\n",(double)Fd);
#endif
#endif

@ @<globals for enumeration@>=
DOUBLE dum1, dum2;

@  @<orthogonalize the basis@>=
#if GIVENS
    givens(lattice,columns,rows,mu,bd,c,N,Fq);  
#else    
    gramschmidt(lattice,columns,rows,mu,bd,c,N,Fq);
#endif    

@ Compute the Fincke-Pohst bounds. See for example H.~Cohen, ``A course
in computational number theory'', page 104.
If we have the quadratic form $Q(x)= x^tR^tRx$, and $r_i$ are the
rows of the matrix $R$  and $r_i'$ are the rows of $R^{-1}$.
Then from $R^{-1}Rx=x$ we obtain $x_i=r_i'Rx$.
With Cauchy-Schwarz this gives
$$
x_i^2 \leq \Vert r_i'\Vert^2(x^tR^tRx)\leq \Vert r_i'\Vert^2C
$$
where $C$ is an upper bound for the solution vector.
We have orthogonalized $B$ such that $B=\tilde{B} \mu$.
And 
$$
\tilde{B}^{-1} = 
\left(\matrix{1/c_1 & 0 & \ldots & 0\cr
0 & 1/c_2 & \ldots & 0 \cr
  & & \ldots & \cr
0 & & \ldots &1/c_n\cr}\right)\cdot \tilde{B}^\top
$$
Therefore, the above $\Vert r_i'\Vert^2$ is in our case
$$r'_{ij} = \sum_l |muinv[i][l]|\cdot|bd[l][j]|/|c[l]|.$$
Then we also compute an upper bound with the $\ell_1$ norm and we take 
the lower of the two bounds.
$$
\vert x_i\vert \leq \Vert r_i'\Vert_1 \Vert Rx\Vert_\infty
$$
@<determine Fincke-Pohst bounds@>=
    fipo_success = 0;
    inverse(mu,muinv,columns);

#ifndef NO_OUTPUT
#if VERBOSE > -1    
    printf("\nFincke-Pohst bounds:\n");@+ fflush(stdout);
#endif    
#endif

    for (i=0;i<columns;i++) {
#if 0   /* Symmetric Fincke-Pohst */
        fipo[i] = 0.0;
        dum1 = 0.0;
        for (j=0;j<rows;j++) {
            dum = 0.0;
            for (l=i;l<columns;l++) dum += muinv[i][l]*bd[l][j]/c[l];
            fipo[i] += dum*dum;
            dum1 += fabs(dum);
        }                
        fipo[i] = SQRT(fipo[i]*Fd);

        dum1 =  fabs(dum1*Fq);
        if (dum1 < fipo[i]) fipo[i] = dum1;

/*        |fipo[i] *= (1.0+EPSILON);|*/

#ifndef NO_OUTPUT
#if VERBOSE > -1    
        printf("%f ",fipo[i]);
#endif
#endif
#else
        /* Asymmetric Fincke-Pohst. */
        fipo[i] = 0.0;
        fipo_u[i] = 0.0;
        fipo_l[i] = 0.0;
        dum1 = 0.0;
        for (j=0;j<rows;j++) {
            for (l=i,dum=0.0;l<columns;l++) dum += muinv[i][l]*bd[l][j]/c[l];
            dum2 = muinv[columns-1][columns-1]*bd[columns-1][j]/c[columns-1];
            fipo[i] += dum*dum;
            dum1 += fabs(dum);

            fipo_u[i] += fabs(dum+dum2);
            fipo_l[i] -= fabs(dum-dum2);
        }                
        fipo[i] = SQRT(fipo[i]*Fd);
        dum1 = fabs(dum1*Fq);
        if (dum1 < fipo[i]) fipo[i] = dum1;

        fipo_u[i] *= Fq;
        fipo_l[i] *= Fq;
#if 0        
        fipo_u[i] -= 1.0;
        fipo_l[i] += 1.0;
        if (fipo[i]<fipo_u[i])  fipo_u[i] = fipo[i];
        if (-fipo[i]>fipo_l[i]) fipo_l[i] = -fipo[i];
#endif        
        fipo_u[i] *= (fipo_u[i]>0)?(1.0+EPSILON):(1.0-EPSILON);
        fipo_l[i] *= (fipo_l[i]>0)?(1.0-EPSILON):(1.0+EPSILON);
#if 1  /* NUMERICAL PROBLEMS!!!! */
        fipo_u[i] = floor(fipo_u[i]);
        fipo_l[i] = ceil(fipo_l[i]);
#endif        
#ifndef NO_OUTPUT
        printf("%03d: %0.3lf\t%0.3lf\n",i, fipo_l[i],fipo_u[i]);
#endif
#endif
    }
#ifndef NO_OUTPUT
#if VERBOSE > -1    
printf("\n");@+ fflush(stdout);
printf("\n");@+ fflush(stdout);
#endif
#endif

@ First, find the index of the first non-zero entry
    in each row.
    Then, detect in each column which rows have its
    first non-zero entry in this column.
@<initialize first-nonzero arrays@>=
    for (l=0;l<rows;l++) {
        for (i=0; i<columns; i++) if (mpz_sgn(get_entry(i,l))!=0) { 
            first_nonzero[l] = i;
            break;
        }           
    }
    j = 0;
    for (l=0;l<columns;l++) {
        firstp[l] = j;
        first_nonzero_in_column[j] = 0;
        j++;
        for (i=0;i<rows;i++) {
            if(first_nonzero[i]==l) {
                first_nonzero_in_column[j] = i; 
                first_nonzero_in_column[firstp[l]]++;
                j++;
            }
        }
    }
    firstp[columns] = j;
    first_nonzero_in_column[j] = 0; 


@ @<globals for enumeration@>=
    long only_zeros_no, only_zeros_success, hoelder_no, hoelder_success;
    long cs_success;
    long N_success;
    
@ @<more initialization@>=
    /*|level = first_nonzero[rows-1];|*/
    level = 0;
    if (level<0) level = 0;
    level_max = level;
    us[level] = 1;
#if 0        
        mpz_set_si(v[level],1);
#else
        v[level] = 1;
#endif
    only_zeros_no = only_zeros_success = 0;
    hoelder_no = hoelder_success = 0;
    cs_success = nosolutions = loops = 0;
    N_success = 0;

@ Increase the loop counter and test if we already can stop after collecting
enough solutions. 
@<increase loop counter@>=
        loops++;
        if ((stop_after_loops>0)&&(stop_after_loops<=loops)) goto afterloop;

#ifndef NO_OUTPUT
#if VERBOSE > -1              
        if (loops%100000000 ==0) {                 /*10000000*/
#if defined(FINCKEPOHST)
            printf("%ld loops, solutions: %ld, fipo: %ld\n",\
            loops,nosolutions, fipo_success);
            /*|
            printf("us: ");
            for (i=columns-1;i>=0;i--) {
                printf("(%d,%d)",i,(int)us[i]);
            }
            printf("\n");
            |*/
#if 0
            /* Write statistics about enumeartion levels */
            for (i=0;i<level_max;i++) {
                printf("%03d: %ld\n", i, nlow[i]);
            }
#endif            
#else            
            printf("%ld loops, solutions: %ld\n",loops,nosolutions);
#endif            
            fflush(stdout);
        }
#endif
#endif

@ The search loop with various pruning tests.
@<the loop of the exhaustive enumeration@>=
    do {  
        @<increase loop counter@>;
        @<compute new |cs|@>;
        if ( (cs[level]<Fd) && (!prune0(fabs(dum),N[level])) )  {   
#if defined(FINCKEPOHST)
#if 0
            if (fabs(us[level])>fipo[level]*(1.0+EPSILON)) {
#else
            if (level!=columns-1 &&
                (us[level]>fipo_u[level] || us[level]<fipo_l[level])
                ) {
#endif                
                fipo_success++;
                goto side_step;
            }            
            
#endif
            compute_w(w,bd,dum,level,rows);

            if (level>0) {
                @<not at a leave@>;
            } else {
                @<at $|level|=0$@>;
            }
        } else {
            cs_success++;
step_back:
/* Up: we go to $|level|+1$. */
            nlow[level]++;
            level++;
            r[level] = level;
            
            if (level_max<level) level_max = level;
side_step:        
/*    
    Side step: the next value in the same level is 
    chosen.
*/
            if (eta[level]==0) {
                delta[level] *= -1;
                if (delta[level]*d[level]>=0) delta[level] += d[level];
            } else {
                delta[level] += (delta[level]*d[level]>=0) ? d[level] : -d[level] ;
            }
#if 0                        
            us[level] = mpz_get_d(v[level]) + delta[level];
#else
            us[level] = v[level] + delta[level];
#endif
        }
    } while (level<columns);     
afterloop:

@ Compute the square of the euclidean length of 
    the projection
@<compute new |cs|@>=
    dum = us[level] + y[level];                      
    cs[level] = cs[level+1] + dum*dum*c[level];
@ We are not at a leave. 
We test, if we can prune the enumeration, otherwise
we decrease the |level|.
@<not at a leave@>=
                if (prune_only_zeros(w[level],level,rows,Fq,first_nonzero_in_column,firstp))
                    goto side_step;

                if ( prune(w[level],cs[level],rows,Fq) ) { 
                    if (eta[level]==1) {
                        goto step_back;
                    }
                    eta[level] = 1;
                    delta[level] *= -1;
                    if (delta[level]*d[level]>=0) delta[level] += d[level];
#if 0                                        
                    us[level] = mpz_get_d(v[level]) + delta[level];
#else
                    us[level] = v[level] + delta[level];
#endif
                } else {
                    level--;
                    eta[level] = 0;
                    delta[level] = 0;
                    if (r[level+1]>r[level]) r[level] = r[level+1];

                    dum = compute_y(mu,us,level,level_max, sigma, r);
                    y[level] = dum;
#if 0                                        
                    mpz_set_d(v[level],ROUND(-dum));
                    us[level] = mpz_get_si(v[level]);
#else
                    v[level] = ROUND(-dum);
                    us[level] = v[level];
#endif
                    d[level] = (ROUND(-dum)>-y[level]) ? -1 : 1; 
                } 
                
@ We arrived at $|level|=0$. We test if we really got a solution and print it.
@<at $|level|=0$@>=
                if (exacttest(w[0],rows,Fq)==1) {          
                    print_solution(w[level],rows,Fq,us,columns); 
#if 0
                    printf("us: ");
                    for (i=columns-1;i>=0;i--) {
                        printf("%d ",(int)us[i]);
                    }
                    printf("after %ld loops\n",loops);
#endif                    
                    if ((stop_after_solutions>0)&&(stop_after_solutions<=nosolutions)) 
                        goto afterloop;
                } 
                goto side_step;
                
@ Additional output at the end of the enumeration.
@<final output@>=
#ifndef NO_OUTPUT
    printf("Prune_cs: %ld\n",cs_success);
    printf("Prune_only_zeros: %ld of %ld\n",only_zeros_success,only_zeros_no);
    printf("Prune_hoelder: %ld of %ld\n",hoelder_success,hoelder_no);
    printf("Prune_N: %ld\n",N_success);
#if defined(FINCKEPOHST)
    printf("Fincke-Pohst: %ld\n",fipo_success); 
#endif    
    printf("Loops: %ld\n",loops);
#endif
    if (( stop_after_solutions<=nosolutions && stop_after_solutions>0 ) ||
         (stop_after_loops<=loops && stop_after_loops>0 )) {
#ifndef NO_OUTPUT
        printf("Stopped after number of solutions: %ld\n",nosolutions);
#endif        
        if ((stop_after_loops<=loops && stop_after_loops>0 )) {
            @<close solution file@>;
            exit(0);
        }
        
    } else {
#ifndef NO_OUTPUT
        printf("Total number of solutions: %ld\n",nosolutions);
#endif
    }    
#ifndef NO_OUTPUT
    printf("\n"); @+ fflush(stdout);
#endif

@ @<free allocated memory for enumeration@>=
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
    for (l=0;l<=columns;l++) free (w[l]);
    free (w);
    free (original_columns);
    for (i=0;i<columns+1;i++) {
       free(sigma[i]);
    }
    free(sigma);
    free(r);
    
#if defined(FINCKEPOHST)
    free (fipo);
    free (fipo_u);
    free (fipo_l);
    for (l=0;l<columns;l++) free (muinv[l]);
    free(muinv);
#endif
#if 0
    free (upb);
    free (lowb);
#endif    
    lllfree(mu,c,N,bd,columns);  

@ @<open solution file@>=
    fp = solfile;
    if (SILENT) fprintf(fp,"SILENT\n");
    fflush(fp);

@ @<close solution file@>=
    if (SILENT) fprintf(fp,"%ld solutions\n",nosolutions);
    fflush(fp);

@ @<global variables@>=
static FILE* fp;

@ @<some vector computations@>=
DOUBLE compute_y(DOUBLE **mu, DOUBLE *us, int level, int level_max, DOUBLE **sigma, int *r) {
    int i;
#if 0
    DOUBLE dum;

    i = level_max;
    dum = 0.0;
    while (i>=level+1) {
        dum += mu[i][level]*us[i];
        i--;
    }
    /*|return dum;|*/
if  (fabs(dum - sigma[level+1][level])>0.000001) {
    printf("diff %0.6lf %0.6lf %0.6lf\n", dum, sigma[level+1][level], fabs(dum - sigma[level+1][level]));
    fflush(stdout);
}
return dum;
#else
    for (i=r[level+1];i>level;i--) sigma[i][level] = sigma[i+1][level] + us[i]*mu[i][level];
    return sigma[level+1][level];
#endif    
}

void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE dum, int level, int rows) {
#if defined(BLAS) 
    dcopy(rows,w[level+1],1,w[level],1);
    daxpy(rows,dum,bd[level],1,w[level],1);
#else
#if defined(BLAS2) 
    cblas_dcopy(rows,w[level+1],1,w[level],1);
    cblas_daxpy(rows,dum,bd[level],1,w[level],1);
#else
    int l;
    
    l = rows-1;
    do {
        w[level][l] = w[level+1][l] + dum*bd[level][l]; 
        l--;
    } while (l>=0);
#endif
#endif
    return;
}

@ The algorithms of Gram-Schmidt and Givens are used
for orthogonalization.
@<orthogonalization@>=
@<Gram-Schmidt algorithm@>;
@<Givens algorithm@>;

@ Orthogonalization with Gram-Schmidt.
@<Gram-Schmidt algorithm@>=
void gramschmidt(COEFF **lattice, int columns, int rows, DOUBLE **mu, 
                 DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq) {
    int i,l,j;
    DOUBLE dum;

    for (i=0; i<columns; i++) {
        for (l=0;l<rows;l++) bd[i][l] = (DOUBLE)mpz_get_d(get_entry(i,l));
        N[i] = 0.0;
        for (j=0;j<i;j++) {
            dum = 0.0;
            for(l=0;l<rows;l++) dum += (DOUBLE)mpz_get_d(get_entry(i,l)) * bd[j][l];
            mu[i][j] = dum / c[j];
            for (l=0;l<rows;l++) bd[i][l] -= mu[i][j]*bd[j][l];
        }
    
        c[i] = scalarproductfp(bd[i],bd[i],rows);
        for(l=0;l<rows;l++) N[i] += fabs(bd[i][l]);
        N[i] /= c[i];
        N[i] *= Fq;
#ifndef NO_OUTPUT
#if VERBOSE > 0     
        printf("%lf ",(double)c[i]); 
#endif
#endif
/*|        N[i] *= (1.0 + EPSILON);|*/
    }
#ifndef NO_OUTPUT
#if VERBOSE > 0
    printf("\n\n"); @+ fflush(stdout);
#endif
#endif
    return;
}

@ Orthogonalization with Givens rotation.
@<Givens algorithm@>=
void givens(COEFF **lattice, int columns, int rows, DOUBLE **mu, 
            DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq) {
    int i,l,j;
    int mm;
    DOUBLE d1,d2;
    DOUBLE gc,gs;
    DOUBLE t;
    

    /* The matrix |b| is copied to |mu|.
       |bd| is set to a $z\times z$ unity matrix.
    */   
    for (i=0; i<columns; i++) {
        for (l=0;l<rows;l++) {
            mu[i][l] = (DOUBLE)mpz_get_d(get_entry(i,l));
        }
    }

    for (i=0;i<rows; i++) {
        for (l=0;l<rows;l++) bd[i][l] = 0.0;
        bd[i][i] = 1.0;
    }

    for (j=1; j<rows; j++) {    /* The Givens rotation */
        mm = (j < columns) ? j : columns;
        for (i=0;i<mm;i++) {
            if (mu[i][j]==0.0) {
                /* Nothing has to be done */
                gc = 1.0;
                gs = 0.0;
            } else {
                /* Stable computation of the
                   rotation coefficients. 
                */
                if (fabs(mu[i][j])>=fabs(mu[i][i])) {
                    t = mu[i][i]/mu[i][j];
                    gs = 1.0 / SQRT(1.0 + t*t);
                    gc = gs * t;
                } else {
                    t = mu[i][j]/mu[i][i];
                    gc = 1.0 / SQRT(1.0 + t*t);
                    gs = gc * t;
                }
                /* Rotation of |mu| */
                for (l=i;l<columns;l++) {        
                    d1 = mu[l][i];
                    d2 = mu[l][j];
                    mu[l][i]  =  gc*d1 + gs*d2;
                    mu[l][j] = -gs*d1 + gc*d2;
                }
                /* Rotation of the matrix $Q^t$ */
                for (l=0;l<rows;l++) {        
                    d1 = bd[i][l];
                    d2 = bd[j][l];
                    bd[i][l]  =  gc*d1 + gs*d2;
                    bd[j][l] = -gs*d1 + gc*d2;
                }
            }
        }
    } 

/* Finally some scaling has to be done, since $Q$ is a orthonormal matrix */
    for (i=0;i<columns;i++) {
        c[i] = mu[i][i]*mu[i][i];
        N[i] = 0.0;
        for (j=0; j<rows; j++) {
            bd[i][j] *= mu[i][i];
            N[i] += fabs(bd[i][j]);
        }
        N[i] /= c[i];
        N[i] *= Fq;

/*|        N[i] *= 1.0 + EPSILON; |*/

        for (j=columns-1; j>=i; j--) mu[j][i] /= mu[i][i];
#ifndef NO_OUTPUT
#if VERBOSE > -1
        printf("%6.3f ",(double)c[i]); 
        if (i>0 && i%15==0) printf("\n");
#endif
#endif
    }
#ifndef NO_OUTPUT
#if VERBOSE > -1
    printf("\n\n"); @+ fflush(stdout);
#endif
#endif
    return;
}
@ The upper triangular matrix $\mu$ is inverted.
@<matrix inversion for Fincke-Pohst@>=
#if defined(FINCKEPOHST)
void inverse(DOUBLE **mu, DOUBLE **muinv, int columns) {
    int i,j,k;
    DOUBLE sum;
    for (j=0;j<columns;j++) 
        for (i=j;i>=0;i--) {
            sum = 0.0;
            for (k=i+1;k<columns;k++) sum += mu[k][i]*muinv[k][j];
            if (i==j) 
                muinv[i][j] = 1.0 - sum;
            else 
                muinv[i][j] = -sum;
        }
    return;
}
#endif

@*There are several pruning methods.
@<pruning subroutines and output@>=
@<exact test@>;
@<simple bound@>;
@<pruning with H\"older@>;
@<prune zeros in row@>;
@<print a solution@>;

@ Exact test of all entries on $|level|=0$.
@<exact test@>=
int exacttest(DOUBLE *v, int rows, DOUBLE Fq) {
    int i;
    i = rows-1;
    do {
        if (fabs(v[i]) > Fq*(1.0 + EPSILON)) { 
            return 0;
        }
        i--;
    } while (i>=0);
    return 1;
}

@ Additional pruning.
@<simple bound@>=
int prune0(DOUBLE li, DOUBLE re) {
    if (li > re*(1+EPSILON)) { 
        N_success++;
        return 1;
    } else {
        return 0;
    }
}

@ Pruning according to H\"olders inequality.
@<pruning with H\"older@>=
int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fq) {
    DOUBLE reseite;
#if !defined(BLAS) 
    int i;
#endif    
    hoelder_no++;
    reseite = (-cs/(1.0+EPSILON))/Fq; /* | * (1-eps) | */

#if defined(BLAS) 
    reseite += dasum(rows,w,1);
#else
#if defined(BLAS2) 
    reseite += cblas_dasum(rows,w,1);
#else
    i = rows-1;
    do {
        reseite += fabs(w[i]); 
        i--;
    } while (i>=0);
#endif    
#endif    
    if (0.0 < reseite) return 0;

    hoelder_success++;
    return 1; 
}

@ Prune if there remain only zeros in a not finished
    row.
@<prune zeros in row@>=
int prune_only_zeros(DOUBLE *w, int level, int rows, DOUBLE Fq, 
                     int *first_nonzero_in_column, int *firstp) {
    int i;
    int f;

    only_zeros_no++;

    if (iszeroone) {
        for (i=0;i<first_nonzero_in_column[firstp[level]];i++) {
            f = first_nonzero_in_column[firstp[level]+1+i];
            if ( fabs(fabs(w[f])-Fq) > 0.5 /* Fq*EPSILON*/ ) {
                only_zeros_success++;
                return 1;
            }
        } 
    } else {
        for (i=0;i<first_nonzero_in_column[firstp[level]];i++) {
            f = first_nonzero_in_column[firstp[level]+1+i];
            if ( fabs(w[f]) > Fq*(1+EPSILON) ) {  
                only_zeros_success++;
                return 1;
            }
        } 
    }
    return 0;
}
@ Output of solutions.
There are difficulties if there are integers bigger than
64 bit integers, even in the 
multiprecision version.
@<print a solution@>=
int print_solution(DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns) {
    int i,j,k;
    int upper;
    int end;
    /* Test again, if the vector is really a solution */
#if 0
    if (fabs(fabs(w[rows-1])-Fq)>0.5*Fq*EPSILON)  {
#else
    if (fabs(w[rows-1])>Fq*(1+EPSILON))  {
#endif  
        return 0; 
    }
    upper = rows-1-free_RHS;
#if 0    
    if (free_RHS && fabs(fabs(w[upper])-Fq)>0.5*Fq*EPSILON) {
#else
    if (free_RHS && fabs(w[upper])>Fq*(1+EPSILON)) {
#endif  
        return 0; 
    }

    if (!SILENT) {
        mpz_set_si(soltest_upfac,1);
#if 0            
        mpz_set_d(soltest_s,ROUND(w[rows-1]));
#else
        mpz_set_si(soltest_s,0);
        for (k=0;k<columns;k++) {
            if (ROUND(us[k])>0) {
                mpz_addmul_ui(soltest_s,get_entry(k,rows-1), ROUND(us[k]));
            } else {
            mpz_submul_ui(soltest_s,get_entry(k,rows-1),-ROUND(us[k]));
            }
        }
#endif            
            
        i = 0;
        if (cut_after_coeff==-1) {
            end=no_original_columns;
#if 0                
            if (nboundvars!=0) { /* Conflicts with |original_columns| */
                end=nboundvars;
            }
#endif                
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
            if (!iszeroone) {
                printf(" ");
                fprintf(fp, " ");
            }
        }  
        if (free_RHS) {
            mpz_set_d(soltest_u,ROUND(w[i]));
            mpz_divexact(soltest_u,soltest_u,max_up);
            mpz_abs(soltest_u,soltest_u);
            printf(" L = ");
            mpz_out_str(NULL,10,soltest_u);
        }
        printf("\n"); @+ fflush(stdout);
        fprintf(fp, "\n"); @+ fflush(fp);
    }

    nosolutions++;
#ifndef NO_OUTPUT
    if (nosolutions%10000==0) {
        printf("%ld\n",nosolutions); @+ fflush(stdout);
    }
#endif

    return 1;
}

@ Output of the polyhedra.
We write the basis vectors as input for the polyhedra programs cdd and lrs.
@<output polyhedra@>=
void basis2poly() {
    return;
}

@ Output of the LP-relaxation.
We write the linear combination of the basis vectors as LP.
@<output LP@>=
void basis2LP(double *low, double *up) {
    return;
}

@*Index.
