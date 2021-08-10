/*4:*/
#line 152 "./diophant.w"

#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h> 
extern long diophant(mpz_t**a_input,mpz_t*b_input,mpz_t*upperbounds_input,
int no_columns,int no_rows,
mpz_t factor_input,mpz_t norm_input,mpz_t scalelastlinefactor,
int silent,int iterate,int iterate_no,
int bkz_beta_input,int bkz_p_input,
long stop_after_sol_input,long stop_after_loops_input,
int free_RHS_input,int*org_col_input,int no_org_col_input,
int cut_after,int nboundedvars,FILE*solfile);

extern void stopProgram();
extern long nosolutions;
#endif

/*:4*//*88:*/
#line 2041 "./diophant.w"

struct constraint{
double val[2];
int parent;
int isSet;
}CONSTRAINT;

/*:88*/
