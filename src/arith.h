#ifndef _ARITH_H
#define _ARITH_H

typedef struct {
    DOUBLE x;
    DOUBLE y;
} hiprec;


extern hiprec twoSum(DOUBLE a, DOUBLE b);
extern void twoSum2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y);
extern void fastTwoSum2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y);
extern void twoSum2i(DOUBLE *a, DOUBLE *b);

extern hiprec split(DOUBLE a);
extern hiprec twoProd(DOUBLE a, DOUBLE b);
extern void twoProd2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y);
extern hiprec twoSquare(DOUBLE a);
extern void twoSquare2(DOUBLE a, DOUBLE *x, DOUBLE *y);

extern DOUBLE hiprec_sqrt(DOUBLE T, DOUBLE t);

extern DOUBLE sumNaive(DOUBLE* p, int n);
extern DOUBLE hiprec_sum2(DOUBLE* p, int n);
extern DOUBLE hiprec_sumK(DOUBLE* p, int n, int K);
extern DOUBLE hiprec_norm_l1(DOUBLE* p, int n);
extern DOUBLE hiprec_normK_l1(DOUBLE* p, int n, int K);

extern DOUBLE dotNaive(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE dotNaiveQP(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE hiprec_dot2(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE hiprec_dotK(DOUBLE* x, DOUBLE* y, int n, int K);
extern DOUBLE hiprec_dot2_row(DOUBLE* x, int dx, DOUBLE* y, int dy, int n);
extern DOUBLE hiprec_normsq_l2(DOUBLE* x, int n);
extern DOUBLE hiprec_norm_l2(DOUBLE* x, int n);
extern DOUBLE hiprec_normK_l2(DOUBLE* x, int n, int K);

#endif