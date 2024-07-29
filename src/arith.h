#ifndef _ARITH_H
#define _ARITH_H

typedef struct {
    DOUBLE hi;
    DOUBLE lo;
} hiprec;

/**
 * Add two doubles and store the result in hi, lo.
 * h, z are helper variables
 */
#define TWOSUM(a, b, hi, lo, h, z) { \
    (h) = (a) + (b);                 \
    (z) = (h) - (a);                 \
    (lo) = ((a) - ((h) - (z))) + ((b) - (z)); \
    (hi) = (h);                      \
}

#define TWOSUM_AVX(a, b, hi, lo, h, z) { \
    (h)  = _mm256_add_pd((a), (b));      \
    (z)  = _mm256_sub_pd((h), (a));      \
    (lo) = _mm256_sub_pd((h), (z));      \
    (lo) = _mm256_sub_pd((a), (lo));     \
    (z)  = _mm256_sub_pd((b), (z));      \
    (lo) = _mm256_add_pd((lo), (z));     \
    (hi) = (h);                          \
}

#define TWOPROD_AVX(a, b, hi, lo) {          \
    (hi)  = _mm256_mul_pd((a), (b));         \
    (lo)  = _mm256_fmsub_pd((a), (b), (hi)); \
}

extern hiprec twoSum(DOUBLE a, DOUBLE b);
extern void twoSum2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo);
extern void fastTwoSum2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo);
extern void twoSum2i(DOUBLE *a, DOUBLE *b);

extern hiprec split(DOUBLE a);
extern hiprec twoProd(DOUBLE a, DOUBLE b);
extern void twoProd2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo);
extern hiprec twoSquare(DOUBLE a);
extern void twoSquare2(DOUBLE a, DOUBLE *x, DOUBLE *y);

extern DOUBLE hiprec_sqrt(DOUBLE T, DOUBLE t);

extern DOUBLE sumNaive(DOUBLE* p, int n);
extern DOUBLE sumNaiveAVX(DOUBLE *p, int n);
extern DOUBLE hiprec_sum2(DOUBLE* p, int n);
extern DOUBLE hiprec_SUM(DOUBLE* p, int n);
extern DOUBLE hiprec_sum_AVX(DOUBLE* p, int n);
extern DOUBLE hiprec_sumK(DOUBLE* p, int n, int K);
extern DOUBLE hiprec_norm_l1(DOUBLE* p, int n);
extern DOUBLE hiprec_normK_l1(DOUBLE* p, int n, int K);
extern DOUBLE hiprec_norm_l1_AVX(DOUBLE* p, int n);

extern DOUBLE dotNaive(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE dotNaiveQP(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE hiprec_dot2(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE hiprec_dot2_AVX(DOUBLE* x, DOUBLE* y, int n);
extern DOUBLE hiprec_dotK(DOUBLE* x, DOUBLE* y, int n, int K);
extern DOUBLE hiprec_dot2_row(DOUBLE* x, int dx, DOUBLE* y, int dy, int n);
extern DOUBLE hiprec_normsq_l2(DOUBLE* x, int n);
extern DOUBLE hiprec_norm_l2(DOUBLE* x, int n);
extern DOUBLE hiprec_normK_l2(DOUBLE* x, int n, int K);

#endif