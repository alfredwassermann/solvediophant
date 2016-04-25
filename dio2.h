#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h>
extern long diophant(mpz_t **a_input, mpz_t *b_input, mpz_t *upperbounds_input,
    int no_columns, int no_rows,
    mpz_t factor_input, mpz_t norm_input, mpz_t scalelastlinefactor,
    int silent, int iterate, int iterate_no,
    int bkz_beta_input, int bkz_p_input,
    long stop_after_sol_input, long stop_after_loops_input,
    int free_RHS_input, int *org_col_input, int no_org_col_input,
    int cut_after, int nboundedvars, FILE* solfile);

extern void stopProgram();
extern long nosolutions;

/* Basic subroutines */
extern void print_num_solutions(long num_solutions);
extern void debug_print(char *m, int l);
extern void print_lattice();
extern long gcd(long n1, long n2);
extern void coeffinit(COEFF *v, int z);
extern int cutlattice();
extern int solutiontest(int position);

#endif
