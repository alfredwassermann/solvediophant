#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h>

#define SQRT sqrt
#define DOUBLE double

/**
 * Definition of the lattice data structures
*/
struct coe {
    mpz_t c;
    int p;
};
#define COEFF struct coe

struct constraint {
    double val[2];
    int parent;
    int isSet;
} CONSTRAINT;


extern long diophant(mpz_t **a_input, mpz_t *b_input, mpz_t *upperbounds_input,
    int no_columns, int no_rows,
    mpz_t factor_input, mpz_t norm_input, mpz_t scalelastlinefactor,
    int silent, int iterate, int iterate_no,
    int bkz_beta_input, int bkz_p_input,
    long stop_after_sol_input, long stop_after_loops_input,
    int free_RHS_input, int *org_col_input, int no_org_col_input,
    int cut_after, int nboundedvars, FILE* solfile);


extern long nosolutions;

extern void stopProgram();
/* Basic subroutines */
extern void print_num_solutions(long num_solutions);
extern void debug_print(char *m, int l);
extern void print_lattice();
extern long gcd(long n1, long n2);
extern void coeffinit(COEFF *v, int z);
extern int cutlattice();
extern int solutiontest(int position);

extern DOUBLE scalarproductlfp (COEFF *v, COEFF *w);
extern DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w , int n);
extern int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N,  DOUBLE ***bs, int s, int z);
extern int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, DOUBLE **bs, int s);
extern double orthogonal_defect(COEFF **lattice, DOUBLE *c, int s, int z);
extern void lll(COEFF **b, int s, int z, DOUBLE quality);
extern DOUBLE iteratedlll(COEFF **b, int s, int z, int no_iterates, DOUBLE quality);
extern DOUBLE bkz(COEFF **b, int s, int z, DOUBLE delta, int beta, int p);
extern DOUBLE enumerate(DOUBLE **mu, DOUBLE *c, long *u, int s, int start_block, int end_block, int p);
extern DOUBLE explicit_enumeration(COEFF **lattice, int columns, int rows);

extern DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max);
extern void compute_w2(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE alpha, int level, int rows);
extern void gramschmidt(COEFF **lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq);
extern void givens(COEFF **lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, DOUBLE Fq);
extern void inverse(DOUBLE **mu, DOUBLE **muinv, int columns);
extern int exacttest(DOUBLE *v, int rows, DOUBLE Fq);
extern int prune0(DOUBLE li, DOUBLE re);
extern int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps);
extern int pruneN(DOUBLE **w, DOUBLE *cs, int t, int rows, int cols, DOUBLE Fq);
extern int prune_only_zeros(DOUBLE **w, int level, int rows, DOUBLE Fq,
                     int *first_nonzero_in_column, int *firstp,
                     DOUBLE **bd, DOUBLE *y, DOUBLE *us, int columns);
extern int prune_snd_nonzero(int columns, int rows,
                     int level, DOUBLE Fq,
                     int *first_nonzero,
                     int *snd_nonzero_in_column, int *sndp,
                     DOUBLE *us,
                     struct constraint *cons);
extern int print_solution(DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns);
extern void shufflelattice();
#endif
