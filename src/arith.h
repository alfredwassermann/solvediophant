/**
Copyright 2024 Alfred Wassermann

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef _ARITH_H
#define _ARITH_H

typedef struct {
    double hi;
    double lo;
} hiprec;

/**
 * Add two doubles and store the result in hi, lo.
 * h, z are helper (double) variables
 */
#define TWOSUM(a, b, hi, lo, h, z) { \
    (h) = (a) + (b);                 \
    (z) = (h) - (a);                 \
    (lo) = ((a) - ((h) - (z))) + ((b) - (z)); \
    (hi) = (h);                      \
}

/**
 * Add two _mm256 words and store the result in hi, lo.
 * h, z are helper (_mm256) variables
 */
#define TWOSUM_AVX(a, b, hi, lo, h, z) { \
    (h)  = _mm256_add_pd((a), (b));      \
    (z)  = _mm256_sub_pd((h), (a));      \
    (lo) = _mm256_sub_pd((h), (z));      \
    (lo) = _mm256_sub_pd((a), (lo));     \
    (z)  = _mm256_sub_pd((b), (z));      \
    (lo) = _mm256_add_pd((lo), (z));     \
    (hi) = (h);                          \
}

/**
 * Add two _mm256 words a and b with a>= b and store the result in hi, lo.
 * h is a helper (_mm256) variable
 */
#define FASTTWOSUM_AVX(a, b, hi, lo, h) { \
    (h)  = _mm256_add_pd((a), (b));       \
    (lo) = _mm256_sub_pd((a), (h));       \
    (lo) = _mm256_add_pd((lo), (b));      \
    (hi) = (h);                           \
}

/**
 * Multiply the _mm256 words a and b using fma.
 * Store the result in hi and lo.
 */
#define TWOPROD_AVX(a, b, hi, lo) {          \
    (hi)  = _mm256_mul_pd((a), (b));         \
    (lo)  = _mm256_fmsub_pd((a), (b), (hi)); \
}

/**
 * Square the _mm256 word a using fma.
 * Store the result in hi and lo.
 */
#define TWOSQUARE_AVX(a, hi, lo) {           \
    (hi)  = _mm256_mul_pd((a), (a));         \
    (lo)  = _mm256_fmsub_pd((a), (a), (hi)); \
}

extern hiprec twoSum(double a, double b);
extern void   twoSum2(double a, double b, double *hi, double *lo);
extern void   fastTwoSum2(double a, double b, double *hi, double *lo);
extern void   twoSum2i(double *a, double *b);

extern hiprec split(double a);
extern hiprec twoProd(double a, double b);
extern void   twoProd2(double a, double b, double *hi, double *lo);
extern hiprec twoSquare(double a);
extern void   twoSquare2(double a, double *x, double *y);

extern double hiprec_sqrt(double T, double t);

extern double hiprec_sum2(double* p, int n);
extern double hiprec_SUM(double* p, int n);
extern double hiprec_sum_AVX(double* p, int n);
extern double hiprec_sumK(double* p, int n, int K);
extern double hiprec_norm_l1(double* p, int n);
extern double hiprec_normK_l1(double* p, int n, int K);
extern double hiprec_norm_l1_AVX(double* p, int n);

extern double dotNaive(double* x, double* y, int n);
extern double dotNaiveQP(double* x, double* y, int n);
extern double hiprec_dot2(double* x, double* y, int n);
extern double hiprec_dot2_AVX(double* x, double* y, int n);
extern double hiprec_dot(double *v, double *w , int n);
extern double double_dot(double *v, double *w , int n);
extern void double_copy(double *to, double *from , int n);
extern void daxpy(double a, double *x, double *y, int n);
extern double double_dot_inc(int n, double *v, int inc_v, double *w , int inc_w);

extern double hiprec_dotK(double* x, double* y, int n, int K);
extern double hiprec_dot2_row(double* x, int dx, double* y, int dy, int n);

extern double hiprec_normsq_l2(double* x, int n);
extern double hiprec_normsq_l2_AVX(double* x, int n);
extern double hiprec_norm_l2(double* x, int n);
extern double hiprec_norm_l2_AVX(double* x, int n);
extern double hiprec_normK_l2(double* x, int n, int K);

extern double hiprec_daxpy_dasum_AVX(double a, double *x, double *y, double *res, int n);
extern double daxpy_dasum_AVX(double a, double *x, double *y, double *res, int n);
extern float saxpy_sasum_AVX(float a, float *x, float *y, float *res, int n);
extern void double_copy_AVX(double* to, double* from, int n);
extern void daxpy_AVX(double a, double *x, double *y, int n);

#if defined(USE_AVX512)
    extern double daxpy_dasum_AVX512(double a, double *x, double *y, double *res, int n);
#endif
#endif