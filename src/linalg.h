#ifndef _LINALG_H
#define _LINALG_H

typedef struct dblexact {
    DOUBLE x;
    DOUBLE y;
} doubleExact;

extern doubleExact twoSum(DOUBLE a, DOUBLE b);
extern void twoSum2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y);
extern void fastTwoSum2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y);
extern void twoSum2i(DOUBLE *a, DOUBLE *b);

extern doubleExact split(DOUBLE a);
extern doubleExact twoProd(DOUBLE a, DOUBLE b);
extern void twoProd2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y);
extern doubleExact twoSquare(DOUBLE a);
extern void twoSquare2(DOUBLE a, DOUBLE *x, DOUBLE *y);
extern DOUBLE accSqrt(DOUBLE T, DOUBLE t);

extern DOUBLE sumNaive(DOUBLE* p, int n);
extern DOUBLE sum2s(DOUBLE* p, int n);
extern DOUBLE sumKvert(DOUBLE* p, int n, int K);
extern DOUBLE sumAbs2s(DOUBLE* p, int n);
extern DOUBLE sumAbsKvert(DOUBLE* p, int n, int K);

extern DOUBLE dotNaive(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE dotNaiveQP(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE dot2(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE norm(DOUBLE* x, int n);

#endif