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
    DOUBLE hi;
    DOUBLE lo;
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

extern hiprec twoSum(DOUBLE a, DOUBLE b);
extern void   twoSum2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo);
extern void   fastTwoSum2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo);
extern void   twoSum2i(DOUBLE *a, DOUBLE *b);

extern hiprec split(DOUBLE a);
extern hiprec twoProd(DOUBLE a, DOUBLE b);
extern void   twoProd2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo);
extern hiprec twoSquare(DOUBLE a);
extern void   twoSquare2(DOUBLE a, DOUBLE *x, DOUBLE *y);

extern DOUBLE hiprec_sqrt(DOUBLE T, DOUBLE t);

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
extern DOUBLE hiprec_normsq_l2_AVX(DOUBLE* x, int n);
extern DOUBLE hiprec_norm_l2(DOUBLE* x, int n);
extern DOUBLE hiprec_norm_l2_AVX(DOUBLE* x, int n);
extern DOUBLE hiprec_normK_l2(DOUBLE* x, int n, int K);

#endif