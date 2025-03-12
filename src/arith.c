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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <x86intrin.h>
// #include <immintrin.h>

#include "const.h"
#include "arith.h"

#if defined(USE_BLAS)
    #define BLAS 1
    #include <cblas-openblas.h>
#elif defined(USE_BLAS_DEV)
    #define BLAS 1
    #include "common.h"
    #include "cblas.h"
#elif defined(USE_BLAS_OLD)
    #define BLAS 1
    #include <cblas.h>
#else
    #define BLAS 0
#endif

/**
* Ogita, Rump, Oishi: Accurate sum and dot product (2005)
*/

/**
 * Add two doubles and return the hiprec result.
 */
hiprec twoSum(double a, double b) {
    hiprec ret;
    double z;

    ret.hi = a + b;
    z = ret.hi - a;
    ret.lo = (a - (ret.hi - z)) + (b - z);

    return ret;
}

/**
 * Add two doubles and return the hiprec result in parameters hi and lo.
 */
void twoSum2(double a, double b, double *hi, double *lo) {
    double z;

    (*hi) = a + b;
    z = (*hi) - a;
    (*lo) = (a - ((*hi) - z)) + (b - z);
}

/**
 * Add two doubles a and b with |a| > |b| and return the hiprec result in parameters hi and lo.
 */
void fastTwoSum2(double a, double b, double *hi, double *lo) {
    // FastTwoSum for |a| > |b|
    (*hi) = a + b;
    (*lo) = (a - (*hi)) + b;
}

/**
 * Add two doubles a and b return hi in a and lo in b.
 */
void twoSum2i(double *a, double *b) {
    /** In-place summation */
    double z, x;

    x = (*a) + (*b);

    z = x - (*a);
    (*b) = ((*a) - (x - z)) + ((*b) - z);
    // FastTwoSum for a > b
    // (*b) = ((*a) - x) + (*b);

    (*a) = x;
}

/**
 * Split a double, see Dekker, and return hiprec
 */
hiprec split(double a) {
    hiprec ret;
    const double factor = 134217729;
    double c;

    c = factor * a;
    ret.hi = (c - (c - a));
    ret.lo = (a - ret.hi);
    return ret;
}

/**
 * Multiply two numbers and return hiprec.
 */
hiprec twoProd(double a, double b) {
    hiprec ret, a_e, b_e;
    ret.hi = a * b;
    a_e = split(a);
    b_e = split(b);
    ret.lo = a_e.lo * b_e.lo - (((ret.hi - a_e.hi * b_e.hi) - a_e.lo * b_e.hi) - a_e.hi * b_e.lo);
    return ret;
}

/**
 * Multiply two numbers and return hiprec result in parameters hi and lo.
 */
void twoProd2(double a, double b, double *hi, double *lo) {
    hiprec a_e, b_e;
    (*hi) = a * b;
    a_e = split(a);
    b_e = split(b);
    (*lo) = a_e.lo * b_e.lo - ((((*hi) - a_e.hi * b_e.hi) - a_e.lo * b_e.hi) - a_e.hi * b_e.lo);
}

/**
 * Multiply two numbers and return hiprec result in parameters hi and lo,
 * use fma()
 */
void twoProd2fma(double a, double b, double *hi, double *lo) {
    (*hi) = a * b;
    (*lo) = fma(a, b, -(*hi));
}

/**
 * Multiply a * a and return hiprec result.
*/
hiprec twoSquare(double a) {
    hiprec ret, a_e;
    ret.hi = a * a;
    a_e = split(a);
    ret.lo = a_e.lo * a_e.lo - ((ret.hi - a_e.hi * a_e.hi) - 2 * a_e.lo * a_e.hi);
    return ret;
}

/**
 * Multiply a * a and return hiprec result in parameters hi and lo.
 */
void twoSquare2(double a, double *hi, double *lo) {
    hiprec a_e;
    (*hi) = a * a;
    a_e = split(a);
    (*lo) = a_e.lo * a_e.lo - (((*hi) - a_e.hi * a_e.hi) - 2 * a_e.lo * a_e.hi);
}

/**
 * Multiply a * a and return hiprec result in parameters hi and lo,
 * use fma().
 */
void twoSquare2fma(double a, double *hi, double *lo) {
    (*hi) = a * a;
    (*lo) = fma(a, a, -(*hi));
}

/**
 *  Accurate square root of hiprec number (T, t).
 */
double hiprec_sqrt(double T, double t) {
    double P, p, H, h;
    double r;

    P = sqrt(T);

    twoSquare2fma(P, &H, &h);
    r = (T - H) - h;
    r = t + r;
    p = r / (2 * P);

    return P + p;
}

/**
 * Sum entries of array p of length n with double precision.
 */
double hiprec_sum2(double* p, int n) {
    double sigma;
    hiprec s;
    int i;

    if (n <= 0) return 0.0;

    s.hi = p[0];
    sigma = 0.0;

    for (i = 1; i < n; i++) {
        s = twoSum(s.hi, p[i]);
        sigma += s.lo;
    }

    return s.hi + sigma;
}

/**
 * Sum entries of array p of length n with double precision.
 */
double hiprec_SUM(double* p, int n) {
    double sigma;
    double h, z;
    hiprec s;
    int i;

    if (n <= 0) return 0.0;

    s.hi = p[0];
    sigma = 0.0;

    for (i = 1; i < n; i++) {
        TWOSUM(s.hi, p[i], s.hi, s.lo, h, z);
        sigma += s.lo;
    }

    return s.hi + sigma;
}

/**
 * hiprec summation using AVX.
 *
 * In general, _mm256_loadu_pd and _mm256_storeu_pd can not be used,
 * since we do not always start to access an array at the beginning.
 * The individual entries of an array are not aligned by 32 bit.
 */
double hiprec_sum_AVX(double* p, int n) {
    hiprec s;
    double sigma_d;
    __m256d vp, lo, h, z;

    if (n <= 0) return 0.0;

    s.hi = 0.0;
    sigma_d = 0.0;

    size_t i = 0;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d sum_hi1 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma1  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi2 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma2  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi3 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma3  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi4 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma4  = {0.0, 0.0, 0.0, 0.0};

    for (i = 0; i < n_16; i += 16) {
        vp = _mm256_loadu_pd(p + i);
        TWOSUM_AVX(sum_hi1, vp, sum_hi1, lo, h, z);
        sigma1 = _mm256_add_pd(sigma1, lo);

        vp = _mm256_loadu_pd(p + i + 4);
        TWOSUM_AVX(sum_hi2, vp, sum_hi2, lo, h, z);
        sigma2 = _mm256_add_pd(sigma2, lo);

        vp = _mm256_loadu_pd(p + i + 8);
        TWOSUM_AVX(sum_hi3, vp, sum_hi3, lo, h, z);
        sigma3 = _mm256_add_pd(sigma3, lo);

        vp = _mm256_loadu_pd(p + i + 12);
        TWOSUM_AVX(sum_hi4, vp, sum_hi4, lo, h, z);
        sigma4 = _mm256_add_pd(sigma4, lo);
    }

    TWOSUM_AVX(sum_hi1, sum_hi2, sum_hi1, lo, h, z);
    sigma1 = _mm256_add_pd(sigma1, sigma2);
    sigma1 = _mm256_add_pd(sigma1, lo);
    TWOSUM_AVX(sum_hi3, sum_hi4, sum_hi3, lo, h, z);
    sigma3 = _mm256_add_pd(sigma3, sigma4);
    sigma3 = _mm256_add_pd(sigma3, lo);
    TWOSUM_AVX(sum_hi1, sum_hi3, sum_hi1, lo, h, z);
    sigma1 = _mm256_add_pd(sigma1, sigma3);
    sigma1 = _mm256_add_pd(sigma1, lo);

    for (i = n_16; i < n_4; i += 4) {
        vp = _mm256_loadu_pd(p + i);
        TWOSUM_AVX(sum_hi1, vp, sum_hi1, lo, h, z);
        sigma1 = _mm256_add_pd(sigma1, lo);
    }

    // Add up the horizontal sum
    s.hi = sum_hi1[0];
    sigma_d = sigma1[0];
    for (i = 1; i < 4; i++) {
        s = twoSum(s.hi, sum_hi1[i]);
        sigma_d += s.lo + sigma1[i];
    }

    // Add the trailing entries
    for (i = n_4; i < n; i++) {
        s = twoSum(s.hi, p[i]);
        sigma_d += s.lo;
    }

    return s.hi + sigma_d;
}

/**
 * Sum entries of array p of length n with K-fold precision.
 */
double hiprec_sumK(double* p, int n, int K) {
    double s, alpha;
    double q[K];
    int i, j, k;

    if (n <= 0) return 0.0;

    K = (n < K) ? n : K;

    for (i = 0; i < K; i++) {
        q[i] = 0.0;
    }

    for (i = 0; i < K - 1; i++) {
        s = p[i];
        for (k = 0; k < i - 1; k++) {
            twoSum2i(&q[k], &s);
        }
        q[i] = s;
    }
    s = q[K - 1];

    for (i = K - 1; i < n; i++) {
        alpha = p[i];
        for (k = 0; k < K - 1; k++) {
            twoSum2i(&q[k], &alpha);
        }
        s += alpha;
    }

    for (j = 0; j < K - 2; j++) {
        alpha = q[j];
        for (k = j + 1; k < K - 1; k++) {
            twoSum2i(&q[k], &alpha);
        }
        s += alpha;
    }

    return s + q[K - 2];
}

/**
 * Sum absolute values of entries of array p of length n with hiprec.
 */
double hiprec_norm_l1(double* p, int n) {
    double sigma;
    hiprec s;
    int i;

    if (n <= 0) return 0.0;

    s.hi = fabs(p[0]);
    sigma = 0.0;

    for (i = 1; i < n; i++) {
        s = twoSum(s.hi, fabs(p[i]));
        sigma += s.lo;
    }

    return s.hi + sigma;
}

/**
 * Sum absolute values of entries of array p of length n with K-fold precision.
 */
double hiprec_normK_l1(double* p, int n, int K) {
    double s, alpha;
    double q[K - 1];
    int i, j, k;

    if (n <= 0 || K < 2) return 0.0;

    K = (n < K) ? n : K;

    for (i = 0; i < K; i++) {
        q[i] = 0.0;
    }

    for (i = 0; i < K - 1; i++) {
        s = fabs(p[i]);
        for (k = 0; k < i - 1; k++) {
            twoSum2i(&q[k], &s);
        }
        q[i] = s;
    }
    s = q[K - 1];

    for (i = K - 1; i < n; i++) {
        alpha = fabs(p[i]);
        for (k = 0; k < K - 1; k++) {
            twoSum2i(&q[k], &alpha);
        }
        s += alpha;
    }

    for (j = 0; j < K - 2; j++) {
        alpha = q[j];
        for (k = j + 1; k < K - 1; k++) {
            twoSum2i(&q[k], &alpha);
        }
        s += alpha;
    }

    return s + q[K - 2];
}

/**
 * Sum absolute values of entries of array p of length n with high precision using AVX.
 */
double hiprec_norm_l1_AVX(double* p, int n) {
    hiprec s;
    double sigma_d;
    __m256d vp, lo, h, z;

    if (n <= 0) return 0.0;

    s.hi = 0.0;
    sigma_d = 0.0;

    size_t i = 0;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d msk = {-0.0, -0.0, -0.0, -0.0}; // Used for fabs()
    __m256d sum_hi1 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma1  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi2 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma2  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi3 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma3  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi4 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma4  = {0.0, 0.0, 0.0, 0.0};

    for (i = 0; i < n_16; i += 16) {
        vp = _mm256_loadu_pd(p + i);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi1, vp, sum_hi1, lo, h, z);
        sigma1 = _mm256_add_pd(sigma1, lo);

        vp = _mm256_loadu_pd(p + i + 4);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi2, vp, sum_hi2, lo, h, z);
        sigma2 = _mm256_add_pd(sigma2, lo);

        vp = _mm256_loadu_pd(p + i + 8);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi3, vp, sum_hi3, lo, h, z);
        sigma3 = _mm256_add_pd(sigma3, lo);

        vp = _mm256_loadu_pd(p + i + 12);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi4, vp, sum_hi4, lo, h, z);
        sigma4 = _mm256_add_pd(sigma4, lo);
    }

    TWOSUM_AVX(sum_hi1, sum_hi2, sum_hi1, lo, h, z);
    sigma1 = _mm256_add_pd(sigma1, sigma2);
    sigma1 = _mm256_add_pd(sigma1, lo);
    TWOSUM_AVX(sum_hi3, sum_hi4, sum_hi3, lo, h, z);
    sigma3 = _mm256_add_pd(sigma3, sigma4);
    sigma3 = _mm256_add_pd(sigma3, lo);
    TWOSUM_AVX(sum_hi1, sum_hi3, sum_hi1, lo, h, z);
    sigma1 = _mm256_add_pd(sigma1, sigma3);
    sigma1 = _mm256_add_pd(sigma1, lo);

    for (i = n_16; i < n_4; i += 4) {
        vp = _mm256_loadu_pd(p + i);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi1, vp, sum_hi1, lo, h, z);
        sigma1 = _mm256_add_pd(sigma1, lo);
    }

    // Add up the horizontal sum
    s.hi = sum_hi1[0];
    sigma_d = sigma1[0];
    for (i = 1; i < 4; i++) {
        s = twoSum(s.hi, sum_hi1[i]);
        sigma_d += s.lo + sigma1[i];
    }

    // Add the trailing entries
    for (i = n_4; i < n; i++) {
        s = twoSum(s.hi, fabs(p[i]));
        sigma_d += s.lo;
    }

    return s.hi + sigma_d;
}

/**
 * Dot product of arrays x and y of length n with high precision.
 */
double hiprec_dot2(double* x, double* y, int n) {
    if (HAS_AVX2) {
        return hiprec_dot2_AVX(x, y, n);
    } else {
        double S, h, s, r, sigma;
        int i;
        if (n <= 0) return 0.0;

        twoProd2fma(x[0], y[0], &S, &sigma);
        for (i = 1; i < n; i++) {
            twoProd2fma(x[i], y[i], &h, &r);
            twoSum2(S, h, &S, &s);
            sigma += (s + r);
        }
        return S + sigma;
    }
}

/**
 * Dot product of arrays x and y of length n with high precision using AVX.
 */
double hiprec_dot2_AVX(double* x, double* y, int n) {
    hiprec s;
    double sigma_d;
    __m256d a, b, hi, lo, q, h, z;

    if (n <= 0) return 0.0;

    s.hi = 0.0;
    sigma_d = 0.0;

    size_t i = 0;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d sum_hi1 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma1  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi2 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma2  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi3 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma3  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi4 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma4  = {0.0, 0.0, 0.0, 0.0};

    for (i = 0; i < n_16; i += 16) {
        a = _mm256_loadu_pd(x + i);
        b = _mm256_loadu_pd(y + i);
        TWOPROD_AVX(a, b, hi, lo);
        TWOSUM_AVX(sum_hi1, hi, sum_hi1, q, h, z);
        lo = _mm256_add_pd(lo, q);
        sigma1 = _mm256_add_pd(sigma1, lo);

        a = _mm256_loadu_pd(x + i + 4);
        b = _mm256_loadu_pd(y + i + 4);
        TWOPROD_AVX(a, b, hi, lo);
        TWOSUM_AVX(sum_hi2, hi, sum_hi2, q, h, z);
        lo = _mm256_add_pd(lo, q);
        sigma2 = _mm256_add_pd(sigma2, lo);

        a = _mm256_loadu_pd(x + i + 8);
        b = _mm256_loadu_pd(y + i + 8);
        TWOPROD_AVX(a, b, hi, lo);
        TWOSUM_AVX(sum_hi3, hi, sum_hi3, q, h, z);
        lo = _mm256_add_pd(lo, q);
        sigma3 = _mm256_add_pd(sigma3, lo);

        a = _mm256_loadu_pd(x + i + 12);
        b = _mm256_loadu_pd(y + i + 12);
        TWOPROD_AVX(a, b, hi, lo);
        TWOSUM_AVX(sum_hi4, hi, sum_hi4, q, h, z);
        lo = _mm256_add_pd(lo, q);
        sigma4 = _mm256_add_pd(sigma4, lo);
    }

    TWOSUM_AVX(sum_hi1, sum_hi2, sum_hi1, lo, h, z);
    lo = _mm256_add_pd(lo, sigma2);
    sigma1 = _mm256_add_pd(sigma1, lo);
    TWOSUM_AVX(sum_hi3, sum_hi4, sum_hi3, lo, h, z);
    lo = _mm256_add_pd(lo, sigma4);
    sigma3 = _mm256_add_pd(sigma3, lo);
    TWOSUM_AVX(sum_hi1, sum_hi3, sum_hi1, lo, h, z);
    lo = _mm256_add_pd(lo, sigma3);
    sigma1 = _mm256_add_pd(sigma1, lo);

    for (i = n_16; i < n_4; i += 4) {
        a = _mm256_loadu_pd(x + i);
        b = _mm256_loadu_pd(y + i);
        TWOPROD_AVX(a, b, hi, lo);
        TWOSUM_AVX(sum_hi1, hi, sum_hi1, q, h, z);
        lo = _mm256_add_pd(lo, q);
        sigma1 = _mm256_add_pd(sigma1, lo);
    }
    // Add up the horizontal sum
    s.hi = sum_hi1[0];
    sigma_d = sigma1[0];
    for (i = 1; i < 4; i++) {
        s = twoSum(s.hi, sum_hi1[i]);
        sigma_d += (s.lo + sigma1[i]);
    }

    // Add the trailing entries
    {
        double hi, lo;
        for (i = n_4; i < n; i++) {
            twoProd2fma(x[i], y[i], &hi, &lo);
            s = twoSum(s.hi, hi);
            sigma_d += (lo + s.lo);
        }
    }

    return s.hi + sigma_d;
}

/**
 * Dot product of arrays x and y of length n with K-fold precision.
*/
double hiprec_dotK(double* x, double* y, int n, int K) {
    double P, H;
    double r[2 * n];
    int i;
    if (n <= 0) return 0.0;

    twoProd2fma(x[0], y[0], &P, &(r[0]));
    for (i = 1; i < n; i++) {
        twoProd2fma(x[i], y[i], &H, &(r[i]));
        twoSum2(P, H, &P, &(r[n + i - 1]));
    }
    r[2 * n - 1] = P;
    return hiprec_sumK(r, 2 * n, K);
}

/**
 * Dot product of arrays x and y of length n and step widths dx and dy with double precision.
 */
double hiprec_dot2_row(double* x, int dx, double* y, int dy, int n) {
    double p, q, s, h, r;
    int i, jx, jy;
    if (n <= 0) return 0.0;

    twoProd2fma(x[0], y[0], &p, &s);
    for (i = 1, jx = dx, jy = dy; i < n; i++, jx += dx, jy += dy) {
        twoProd2fma(x[jx], y[jy], &h, &r);
        twoSum2(p, h, &p, &q);
        s += (q + r);
    }
    return p + s;
}

/**
 * @private
 * Square of ||x||_2, i.e.
 * dot product of array x with itself of length n with high precision.
 * Used in hiprec_normsq_l2_AVX and hiprec_norm_l2_AVX.
 * @returns hiprec
 *
 */
hiprec hiprec_normsq_l2_kernel(double* x, int n) {
    double P, p, H, h;
    double d;
    hiprec s;
    int i;

    s.hi = 0.0;
    s.lo = 0.0;

    if (n <= 0) return s;

    for (i = 0; i < n; i++) {
        twoSquare2fma(x[i], &P, &p);
        fastTwoSum2(s.hi, P, &H, &h);
        d = h + (s.lo + p);
        fastTwoSum2(H, d, &(s.hi), &(s.lo));
    }
    return s;
}

/**
 * Square of ||x||_2, i.e.
 * dot product of array x with itself of length n with high precision.
 */
double hiprec_normsq_l2(double* x, int n) {
    if (HAS_AVX2) {
        return hiprec_normsq_l2_AVX(x, n);
    } else {
        hiprec s = hiprec_normsq_l2_kernel(x, n);
        return s.hi + s.lo;
    }
}

/**
 * ||x||_2, i.e. square root of dot product of array x with itself of length n with high precision.
 * "Fast and accurate computation of the Euclidean norm of a vector"
 * Siegfried M. Rump (2007)
 */
double hiprec_norm_l2(double* x, int n) {
    if (HAS_AVX2) {
        return hiprec_norm_l2_AVX(x, n);
    } else {
        hiprec s = hiprec_normsq_l2_kernel(x, n);
        return hiprec_sqrt(s.hi, s.lo);
    }
}

/**
 * @private
 * Square of ||x||_2, i.e.
 * dot product of array x with itself of length n with high precision using AVX.
 * Used in hiprec_normsq_l2_AVX and hiprec_norm_l2_AVX.
 * @returns hiprec
 *
 */
hiprec hiprec_normsq_l2_AVX_kernel(double* x, int n) {
    hiprec s;
    double sigma_d;
    __m256d a, hi, lo, q, h, z;

    s.hi = 0.0;
    s.lo = 0.0;
    sigma_d = 0.0;

    if (n <= 0) return s;

    size_t i = 0;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d sum_hi1 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma1  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi2 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma2  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi3 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma3  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi4 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma4  = {0.0, 0.0, 0.0, 0.0};

    for (i = 0; i < n_16; i += 16) {
        a = _mm256_loadu_pd(x + i);
        TWOSQUARE_AVX(a, hi, lo);
        // TWOSUM_AVX(sum_hi1, hi, sum_hi1, q, h, z);
        FASTTWOSUM_AVX(sum_hi1, hi, sum_hi1, q, h);
        lo = _mm256_add_pd(lo, q);
        sigma1 = _mm256_add_pd(sigma1, lo);

        a = _mm256_loadu_pd(x + i + 4);
        TWOSQUARE_AVX(a, hi, lo);
        // TWOSUM_AVX(sum_hi2, hi, sum_hi2, q, h, z);
        FASTTWOSUM_AVX(sum_hi2, hi, sum_hi2, q, h);
        lo = _mm256_add_pd(lo, q);
        sigma2 = _mm256_add_pd(sigma2, lo);

        a = _mm256_loadu_pd(x + i + 8);
        TWOSQUARE_AVX(a, hi, lo);
        // TWOSUM_AVX(sum_hi3, hi, sum_hi3, q, h, z);
        FASTTWOSUM_AVX(sum_hi3, hi, sum_hi3, q, h);
        lo = _mm256_add_pd(lo, q);
        sigma3 = _mm256_add_pd(sigma3, lo);

        a = _mm256_loadu_pd(x + i + 12);
        TWOSQUARE_AVX(a, hi, lo);
        // TWOSUM_AVX(sum_hi4, hi, sum_hi4, q, h, z);
        FASTTWOSUM_AVX(sum_hi4, hi, sum_hi4, q, h);
        lo = _mm256_add_pd(lo, q);
        sigma4 = _mm256_add_pd(sigma4, lo);
    }

    TWOSUM_AVX(sum_hi1, sum_hi2, sum_hi1, lo, h, z);
    lo = _mm256_add_pd(lo, sigma2);
    sigma1 = _mm256_add_pd(sigma1, lo);
    TWOSUM_AVX(sum_hi3, sum_hi4, sum_hi3, lo, h, z);
    lo = _mm256_add_pd(lo, sigma4);
    sigma3 = _mm256_add_pd(sigma3, lo);
    TWOSUM_AVX(sum_hi1, sum_hi3, sum_hi1, lo, h, z);
    lo = _mm256_add_pd(lo, sigma3);
    sigma1 = _mm256_add_pd(sigma1, lo);

    for (i = n_16; i < n_4; i += 4) {
        a = _mm256_loadu_pd(x + i);
        TWOSQUARE_AVX(a, hi, lo);
        // TWOSUM_AVX(sum_hi1, hi, sum_hi1, q, h, z);
        FASTTWOSUM_AVX(sum_hi1, hi, sum_hi1, q, h);
        lo = _mm256_add_pd(lo, q);
        sigma1 = _mm256_add_pd(sigma1, lo);
    }

    // Add up the horizontal sum
    s.hi = sum_hi1[0];
    sigma_d = sigma1[0];
    for (i = 1; i < 4; i++) {
        s = twoSum(s.hi, sum_hi1[i]);
        sigma_d += s.lo + sigma1[i];
    }

    // Add the trailing entries
    {
        double hi, lo;
        for (i = n_4; i < n; i++) {
            twoSquare2fma(x[i], &hi, &lo);
            s = twoSum(s.hi, hi);
            sigma_d += (lo + s.lo);
        }
    }
    s.lo = sigma_d;
    return s;
}

/**
 * Square of ||x||_2, i.e.
 * dot product of array x with itself of length n with high precision using AVX.
 */
double hiprec_normsq_l2_AVX(double* x, int n) {
    hiprec s = hiprec_normsq_l2_AVX_kernel(x, n);
    return s.hi + s.lo;
}

/**
 * ||x||_2, i.e. square root of dot product of array x with itself of length n with high precision
 * using AVX.
 */
double hiprec_norm_l2_AVX(double* x, int n) {
    if (n <= 0) {
        return 0.0;
    }
    hiprec s = hiprec_normsq_l2_AVX_kernel(x, n);
    return hiprec_sqrt(s.hi, s.lo);
}

/**
 * ||x||_2, i.e. square root of dot product of array x with itself of length n with K-fold precision.
 */
double hiprec_normK_l2(double* x, int n, int K) {
    double S, P;
    double r[2 * n];
    int i;

    if (n <= 0) return 0.0;

    S = 0.0;
    for (i = 0; i < n; i++) {
        twoSquare2fma(x[i], &P, &(r[i]));
        twoSum2(S, P, &S, &(r[n + i - 1]));
    }
    r[2 * n - 1] = S;

    return hiprec_sqrt(hiprec_sumK(r, 2 * n, K), 0.0);
}

/**
 * Combine daxpy and dasum in one loop, hiprec version.
 * Unused.
 */
double hiprec_daxpy_dasum_AVX(double a, double *x, double *y, double *res, int n) {
    hiprec s;
    double sigma_d;
    __m256d va, vx, vy, vp, lo, h, z;

    if (n <= 0) return 0.0;

    size_t i = 0;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    s.hi = 0.0;
    sigma_d = 0.0;

    __m256d msk = {-0.0, -0.0, -0.0, -0.0}; // Used for fabs()

    __m256d sum_hi1 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma1  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi2 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma2  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi3 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma3  = {0.0, 0.0, 0.0, 0.0};
    __m256d sum_hi4 = {0.0, 0.0, 0.0, 0.0};
    __m256d sigma4  = {0.0, 0.0, 0.0, 0.0};

    va = _mm256_broadcast_sd(&a);

    for (i = 0; i < n_16; i += 16) {
        vx = _mm256_loadu_pd(x + i);
        vy = _mm256_loadu_pd(y + i);
        vp = _mm256_fmadd_pd(va, vx, vy);
        _mm256_storeu_pd(res + i, vp);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi1, vp, sum_hi1, lo, h, z);
        sigma1 = _mm256_add_pd(sigma1, lo);

        vx = _mm256_loadu_pd(x + i + 4);
        vy = _mm256_loadu_pd(y + i + 4);
        vp = _mm256_fmadd_pd(va, vx, vy);
        _mm256_storeu_pd(res + i + 4, vp);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi2, vp, sum_hi2, lo, h, z);
        sigma2 = _mm256_add_pd(sigma2, lo);

        vx = _mm256_loadu_pd(x + i + 8);
        vy = _mm256_loadu_pd(y + i + 8);
        vp = _mm256_fmadd_pd(va, vx, vy);
        _mm256_storeu_pd(res + i + 8, vp);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi3, vp, sum_hi3, lo, h, z);
        sigma3 = _mm256_add_pd(sigma3, lo);

        vx = _mm256_loadu_pd(x + i + 12);
        vy = _mm256_loadu_pd(y + i + 12);
        vp = _mm256_fmadd_pd(va, vx, vy);
        _mm256_storeu_pd(res + i + 12, vp);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi4, vp, sum_hi4, lo, h, z);
        sigma4 = _mm256_add_pd(sigma4, lo);
    }

    TWOSUM_AVX(sum_hi1, sum_hi2, sum_hi1, lo, h, z);
    sigma1 = _mm256_add_pd(sigma1, sigma2);
    sigma1 = _mm256_add_pd(sigma1, lo);
    TWOSUM_AVX(sum_hi3, sum_hi4, sum_hi3, lo, h, z);
    sigma3 = _mm256_add_pd(sigma3, sigma4);
    sigma3 = _mm256_add_pd(sigma3, lo);
    TWOSUM_AVX(sum_hi1, sum_hi3, sum_hi1, lo, h, z);
    sigma1 = _mm256_add_pd(sigma1, sigma3);
    sigma1 = _mm256_add_pd(sigma1, lo);

    for (i = n_16; i < n_4; i += 4) {
        vx = _mm256_loadu_pd(x + i);
        vy = _mm256_loadu_pd(y + i);
        vp = _mm256_fmadd_pd(va, vx, vy);
        _mm256_storeu_pd(res + i, vp);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        TWOSUM_AVX(sum_hi1, vp, sum_hi1, lo, h, z);
        sigma1 = _mm256_add_pd(sigma1, lo);
    }

    // Add up the horizontal sum
    s.hi = sum_hi1[0];
    sigma_d = sigma1[0];
    for (i = 1; i < 4; i++) {
        s = twoSum(s.hi, sum_hi1[i]);
        sigma_d += s.lo + sigma1[i];
    }

    // Handle the trailing entries
    for (i = n_4; i < n; i++) {
        res[i] = fma(a, x[i], y[i]);
        s = twoSum(s.hi, fabs(res[i]));
        sigma_d += s.lo;
    }

    return s.hi + sigma_d;
}

/**
 * Combine daxpy and dasum in one loop, feeding 4 pipes.
 * Seems to be the fastest one. Use AVX2.
 */
double daxpy_dasum_AVX(double a, double *x, double *y, double *res, int n) {
    size_t i = 0;
    double s = 0.0;
    __m256d va, vp1, vp2, vp3, vp4;

    if (n <= 0) return 0.0;

    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d msk = {-0.0, -0.0, -0.0, -0.0}; // Used for fabs()

    __m256d sum1 = {0.0, 0.0, 0.0, 0.0};
    __m256d sum2 = {0.0, 0.0, 0.0, 0.0};
    __m256d sum3 = {0.0, 0.0, 0.0, 0.0};
    __m256d sum4 = {0.0, 0.0, 0.0, 0.0};

    va = _mm256_broadcast_sd(&a);

    /*
    if (((long)x & 31) != 0 || ((long)y & 31) != 0) {
        // Address is not 32 bit aligned
        fprintf(stderr, "Unaligned array in daxpy_dasum_AVX\n");
        if (( (long)x & 31) != 0) {
            fprintf(stderr, "x is unaligned array %ld\n", (size_t)x);
        }
        if (( (long)y & 31) != 0) {
            fprintf(stderr, "y is unaligned array\n");
        }
        exit(EXIT_MEMORY);
    }
    */

    for (i = 0; i < n_16; i += 16) {
        vp1 = _mm256_fmadd_pd(va, _mm256_load_pd(x),      _mm256_load_pd(y));
        vp2 = _mm256_fmadd_pd(va, _mm256_load_pd(x + 4),  _mm256_load_pd(y + 4));
        vp3 = _mm256_fmadd_pd(va, _mm256_load_pd(x + 8),  _mm256_load_pd(y + 8));
        vp4 = _mm256_fmadd_pd(va, _mm256_load_pd(x + 12), _mm256_load_pd(y + 12));

        sum1 = _mm256_add_pd(sum1, _mm256_andnot_pd(msk, vp1));
        sum2 = _mm256_add_pd(sum2, _mm256_andnot_pd(msk, vp2));
        sum3 = _mm256_add_pd(sum3, _mm256_andnot_pd(msk, vp3));
        sum4 = _mm256_add_pd(sum4, _mm256_andnot_pd(msk, vp4));

        _mm256_store_pd(res, vp1);
        _mm256_store_pd(res + 4, vp2);
        _mm256_store_pd(res + 8, vp3);
        _mm256_store_pd(res + 12, vp4);

        res += 16;
        x += 16;
        y += 16;
    }

    sum1 = _mm256_add_pd(sum1, sum2);
    sum3 = _mm256_add_pd(sum3, sum4);
    sum1 = _mm256_add_pd(sum1, sum3);

    for (i = n_16; i < n_4; i += 4) {
        vp1 = _mm256_fmadd_pd(va, _mm256_load_pd(x), _mm256_load_pd(y));
        _mm256_store_pd(res, vp1);
        sum1 = _mm256_add_pd(sum1, _mm256_andnot_pd(msk, vp1));

        res += 4;
        x += 4;
        y += 4;
    }

    // Add up the horizontal sum
    s = sum1[0] + sum1[1] + sum1[2] + sum1[3];

    // Handle the trailing entries
    for (i = n_4; i < n; i++) {
        *res = fma(a, *x, *y);
        s += fabs(*res);
        res++;
        x++;
        y++;
    }

    return s;
}


/**
 * Combine saxpy and sasum in one loop, feeding 4 pipes.
 * Seems to be the fastest one. Use AVX2.
 */
float saxpy_sasum_AVX(float a, float *x, float *y, float *res, int n) {
    size_t i = 0;
    float s = 0.0;
    __m256 va, vp1, vp2, vp3, vp4;

    if (n <= 0) return 0.0;

    size_t n_32 = n & -32;
    size_t n_8 = n & -8;

    __m256 msk = {-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0};  // Used for fabs()
    __m256 sum1 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    __m256 sum2 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    __m256 sum3 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    __m256 sum4 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    va = _mm256_broadcast_ss(&a);

    for (i = 0; i < n_32; i += 32) {
        vp1 = _mm256_fmadd_ps(va, _mm256_load_ps(x),      _mm256_load_ps(y));
        vp2 = _mm256_fmadd_ps(va, _mm256_load_ps(x + 8),  _mm256_load_ps(y + 8));
        vp3 = _mm256_fmadd_ps(va, _mm256_load_ps(x + 16), _mm256_load_ps(y + 16));
        vp4 = _mm256_fmadd_ps(va, _mm256_load_ps(x + 24), _mm256_load_ps(y + 24));

        sum1 = _mm256_add_ps(sum1, _mm256_andnot_ps(msk, vp1));
        sum2 = _mm256_add_ps(sum2, _mm256_andnot_ps(msk, vp2));
        sum3 = _mm256_add_ps(sum3, _mm256_andnot_ps(msk, vp3));
        sum4 = _mm256_add_ps(sum4, _mm256_andnot_ps(msk, vp4));

        _mm256_store_ps(res,      vp1);
        _mm256_store_ps(res + 8,  vp2);
        _mm256_store_ps(res + 16, vp3);
        _mm256_store_ps(res + 24, vp4);

        res += 32;
        x += 32;
        y += 32;
    }

    sum1 = _mm256_add_ps(sum1, sum2);
    sum3 = _mm256_add_ps(sum3, sum4);
    sum1 = _mm256_add_ps(sum1, sum3);

    for (i = n_32; i < n_8; i += 8) {
        vp1 = _mm256_fmadd_ps(va, _mm256_load_ps(x), _mm256_load_ps(y));
        sum1 = _mm256_add_ps(sum1, _mm256_andnot_ps(msk, vp1));
        _mm256_store_ps(res, vp1);

        res += 8;
        x += 8;
        y += 8;
    }

    // Add up the horizontal sum
    s = sum1[0] + sum1[1] + sum1[2] + sum1[3] + sum1[4] + sum1[5] + sum1[6] + sum1[7];
    
    // Handle the trailing entries
    for (i = n_8; i < n; i++) {
        *res = fma(a, *x, *y);
        s += fabs(*res);
        res++;
        x++;
        y++;
    }

    return s;
}

void daxpy_AVX(double a, double *x, double *y, int n) {
    if (n <= 0) return;

    size_t i = 0;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d va = _mm256_broadcast_sd(&a);

    for (i = 0; i < n_16; i += 16) {
        _mm256_storeu_pd(y,      _mm256_fmadd_pd(va, _mm256_loadu_pd(x),      _mm256_loadu_pd(y)));
        _mm256_storeu_pd(y + 4,  _mm256_fmadd_pd(va, _mm256_loadu_pd(x + 4),  _mm256_loadu_pd(y + 4)));
        _mm256_storeu_pd(y + 8,  _mm256_fmadd_pd(va, _mm256_loadu_pd(x + 8),  _mm256_loadu_pd(y + 8)));
        _mm256_storeu_pd(y + 12, _mm256_fmadd_pd(va, _mm256_loadu_pd(x + 12), _mm256_loadu_pd(y + 12)));
        x += 16;
        y += 16;
    }

    for (i = n_16; i < n_4; i += 4) {
        _mm256_storeu_pd(y,      _mm256_fmadd_pd(va, _mm256_loadu_pd(x),      _mm256_loadu_pd(y)));
        x += 4;
        y += 4;
    }

    // Handle the trailing entries
    for (i = n_4; i < n; i++) {
        *y = fma(a, *x, *y);
        x++;
        y++;
    }
}

double double_dot_AVX(double* x, double* y, int n) {
    double s = 0.0;
    __m256d a, b;

    if (n <= 0) return 0.0;

    size_t i = 0;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d sum1 = {0.0, 0.0, 0.0, 0.0};
    __m256d sum2 = {0.0, 0.0, 0.0, 0.0};
    __m256d sum3 = {0.0, 0.0, 0.0, 0.0};
    __m256d sum4 = {0.0, 0.0, 0.0, 0.0};

    for (i = 0; i < n_16; i += 16) {
        a = _mm256_loadu_pd(x);
        b = _mm256_loadu_pd(y);
        sum1 = _mm256_fmadd_pd(a, b, sum1);

        a = _mm256_loadu_pd(x + 4);
        b = _mm256_loadu_pd(y + 4);
        sum2 = _mm256_fmadd_pd(a, b, sum2);

        a = _mm256_loadu_pd(x + 8);
        b = _mm256_loadu_pd(y + 8);
        sum3 = _mm256_fmadd_pd(a, b, sum3);

        a = _mm256_loadu_pd(x + 12);
        b = _mm256_loadu_pd(y + 12);
        sum4 = _mm256_fmadd_pd(a, b, sum4);

        x += 16;
        y += 16;
    }

    sum1 = _mm256_add_pd(sum1, sum2);
    sum3 = _mm256_add_pd(sum3, sum4);
    sum1 = _mm256_add_pd(sum1, sum3);

    for (i = n_16; i < n_4; i += 4) {
        a = _mm256_loadu_pd(x);
        b = _mm256_loadu_pd(y);
        sum1 = _mm256_fmadd_pd(a, b, sum1);
        x += 4;
        y += 4;
    }

    // Add up the horizontal sum
    s = sum1[0] + sum1[1] + sum1[2] + sum1[3];

    // Add the trailing entries
    for (i = n_4; i < n; i++) {
        s = fma(*x, *y, s);
        x++;
        y++;
    }

    return s;
}

void double_copy_AVX(double* to, double* from, int n) {
    if (n <= 0) return;

    size_t i;
    size_t n_16 = n & -16;
    size_t n_4 = n & -4;

    __m256d f1, f2, f3, f4;

    for (i = 0; i < n_16; i += 16) {
        f1 = _mm256_loadu_pd(from + i);
        f2 = _mm256_loadu_pd(from + i + 4);
        f3 = _mm256_loadu_pd(from + i + 8);
        f4 = _mm256_loadu_pd(from + i + 12);

        _mm256_storeu_pd(to + i, f1);
        _mm256_storeu_pd(to + i + 4, f2);
        _mm256_storeu_pd(to + i + 8, f3);
        _mm256_storeu_pd(to + i + 12, f4);
    }

    for (i = n_16; i < n_4; i += 4) {
        f1 = _mm256_loadu_pd(from + i);
        _mm256_storeu_pd(to + i, f1);
    }

    // Do the trailing entries
    {
        for (i = n_4; i < n; i++) {
            to[i] = from[i];
        }
    }
}

/**
 * Combine daxpy and dasum in one loop, feeding 4 pipes.
 * Seems to be the fastest one. Use AVX2.
 */
#if defined(USE_AVX512)
double daxpy_dasum_AVX512(double a, double *x, double *y, double *res, int n) {
    size_t i = 0;
    double s = 0.0;
    __m512d va, vp1, vp2, vp3, vp4;

    if (n <= 0) return 0.0;

    size_t n_32 = n & -32;

    if (n_32 >= 32) {
        __m512d msk = {-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0}; // Used for fabs()

        __m512d sum1 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        __m512d sum2 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        __m512d sum3 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        __m512d sum4 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        // va = _mm512_broadcast_pd(&a);
        va = _mm512_set1_pd(a);

        if (((long)x & 31) != 0 || ((long)y & 31) != 0) {
            // Address is not 32 bit aligned
            fprintf(stderr, "Unaligned array in daxpy_dasum_AVX512\n");
            if (( (long)x & 31) != 0) {
                fprintf(stderr, "x is unaligned array %ld\n", (size_t)x);
            }
            if (( (long)y & 31) != 0) {
                fprintf(stderr, "y is unaligned array\n");
            }
            exit(EXIT_MEMORY);
        }

        for (i = 0; i < n_32; i += 32) {
            vp1 = _mm512_fmadd_pd(va, _mm512_load_pd(x + i), _mm512_load_pd(y + i));
            _mm512_store_pd(res + i, vp1);
            sum1 = _mm512_add_pd(sum1, _mm512_andnot_pd(msk, vp1));     // vp = _mm512_andnot_pd(msk, vp); // fabs(vp)

            vp2 = _mm512_fmadd_pd(va, _mm512_load_pd(x + i + 4), _mm512_load_pd(y + i + 4));
            _mm512_store_pd(res + i + 4, vp2);
            sum2 = _mm512_add_pd(sum2, _mm512_andnot_pd(msk, vp2));

            vp3 = _mm512_fmadd_pd(va, _mm512_load_pd(x + i + 8), _mm512_load_pd(y + i + 8));
            _mm512_store_pd(res + i + 8, vp3);
            sum3= _mm512_add_pd(sum3, _mm512_andnot_pd(msk, vp3));

            vp4 = _mm512_fmadd_pd(va, _mm512_load_pd(x + i + 12), _mm512_load_pd(y + i + 12));
            _mm512_store_pd(res + i + 12, vp4);
            sum4 = _mm512_add_pd(sum4, _mm512_andnot_pd(msk, vp4));
        }

        sum1 = _mm512_add_pd(sum1, sum2);
        sum3 = _mm512_add_pd(sum3, sum4);
        sum1 = _mm512_add_pd(sum1, sum3);

        // Add up the horizontal sum
        s = sum1[0];
        for (i = 1; i < 8; i++) {
            s += sum1[i];
        }
    }

    // Handle the trailing entries
    for (i = n_32; i < n; i++) {
        res[i] = fma(a, x[i], y[i]);
        s += fabs(res[i]);
    }

    return s;
}
#endif

double hiprec_dot(double *v, double *w , int n) {
    if (HAS_AVX2) {
        return hiprec_dot2_AVX(v, w, n);
    } else {
        return hiprec_dot2(v, w, n);
    }
}

double double_dot(double *v, double *w , int n) {
    if (HAS_AVX2) {
        return double_dot_AVX(v, w, n);
    } else {
        #if BLAS
            return cblas_ddot(n, v, 1, w, 1);
        #else
            double r;
            int i;
            r = 0.0;
            for (i = n - 1; i >= 0; i--) r += v[i] * w[i];
            return r;
        #endif
    }
}

double double_dot_inc(int n, double *v, int inc_v, double *w , int inc_w) {
    int i = 0;
    int iv = 0, iw = 0;
    double s = 0.0;

    if (n <= 0) return s;

    double tmp1 = 0.0;
    double tmp2 = 0.0;

    int n1 = n & -4;

    while (i < n1) {
        double m1 = w[iw]             * v[iv] ;
        double m2 = w[iw + inc_w]     * v[iv + inc_v];
        double m3 = w[iw + 2 * inc_w] * v[iv + 2 * inc_v];
        double m4 = w[iw + 3 * inc_w] * v[iv + 3 * inc_v];

        iv += inc_v * 4;
        iw += inc_w * 4;

        tmp1 += m1 + m3;
        tmp2 += m2 + m4;

        i += 4;
    }

    while(i < n) {
        tmp1 += w[iw] * v[iv] ;
        iv  += inc_v;
        iw  += inc_w;
        i++ ;
    }
    s = tmp1 + tmp2;
    return s;
}

void double_copy(double *to, double *from , int n) {
    if (HAS_AVX2) {
        double_copy_AVX(to, from, n);
    } else {
        #if BLAS
            cblas_dcopy(n, from, 1, to, 1);
        #else
            for (int i = 0; i < n; ++i) {
                to[i] = from[i];
            }
        #endif
    }
}

void daxpy(double a, double *x, double *y, int n) {
    if (HAS_AVX2) {
        daxpy_AVX(a, x, y, n);
    } else {
        #if BLAS
            cblas_daxpy(n, a, x, 1, y, 1);
        #else
            for (int j = 0; j < n; ++j) {
                x[j] += a * y[j];
            }
        #endif
    }
}

/* -------------  Experiments with number of pipes --------------------- */

/**
 * Combine daxpy and dasum in one loop, simple loop.
 * Slower than 4 pipes.
 */
double daxpy_dasum_AVX1(double a, double *x, double *y, double *res, int n) {
    size_t i = 0;
    double s = 0.0;
    __m256d va, vx, vy, vp;

    if (n <= 0) return 0.0;

    size_t n_4 = n & -4;

    __m256d msk = {-0.0, -0.0, -0.0, -0.0}; // Used for fabs()
    __m256d sum1 = {0.0, 0.0, 0.0, 0.0};

    va = _mm256_broadcast_sd(&a);

    if (((long)x & 31) != 0 || ((long)y & 31) != 0) {
        // Address is not 32 bit aligned
        fprintf(stderr, "Unaligned array in daxpy_dasum_AVX1\n");
        exit(EXIT_MEMORY);
    }

    for (i = 0; i < n_4; i += 4) {
        vx = _mm256_load_pd(x + i);
        vy = _mm256_load_pd(y + i);
        vp = _mm256_fmadd_pd(va, vx, vy);
        _mm256_store_pd(res + i, vp);
        vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
        sum1 = _mm256_add_pd(sum1, vp);
    }

    // Add up the horizontal sum
    s = sum1[0];
    for (i = 1; i < 4; i++) {
        s += sum1[i];
    }

    // Handle the trailing entries
    for (i = n_4; i < n; i++) {
        res[i] = fma(a, x[i], y[i]);
        s += fabs(res[i]);
    }

    return s;
}

/**
 * Combine daxpy and dasum in one loop, feeding 8 pipes.
 * Slower than 4 pipes.
 */
double daxpy_dasum_AVX8(double a, double *x, double *y, double *res, int n) {
    size_t i = 0;
    double s = 0.0;
    __m256d va, vx, vy, vp;

    if (n <= 0) return 0.0;

    size_t n_32 = n & -32;

    if (n_32 >= 32) {
        __m256d msk = {-0.0, -0.0, -0.0, -0.0}; // Used for fabs()

        __m256d sum1 = {0.0, 0.0, 0.0, 0.0};
        __m256d sum2 = {0.0, 0.0, 0.0, 0.0};
        __m256d sum3 = {0.0, 0.0, 0.0, 0.0};
        __m256d sum4 = {0.0, 0.0, 0.0, 0.0};
        __m256d sum5 = {0.0, 0.0, 0.0, 0.0};
        __m256d sum6 = {0.0, 0.0, 0.0, 0.0};
        __m256d sum7 = {0.0, 0.0, 0.0, 0.0};
        __m256d sum8 = {0.0, 0.0, 0.0, 0.0};

        va = _mm256_broadcast_sd(&a);
        if (((long)x & 31) != 0 || ((long)y & 31) != 0) {
            // Address is not 32 bit aligned
            fprintf(stderr, "Unaligned array in daxpy_dasum_AVX8\n");
            exit(EXIT_MEMORY);
        }

        for (i = 0; i < n_32; i += 32) {
            vx = _mm256_load_pd(x + i);
            vy = _mm256_load_pd(y + i);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum1 = _mm256_add_pd(sum1, vp);

            vx = _mm256_load_pd(x + i + 4);
            vy = _mm256_load_pd(y + i + 4);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i + 4, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum2 = _mm256_add_pd(sum2, vp);

            vx = _mm256_load_pd(x + i + 8);
            vy = _mm256_load_pd(y + i + 8);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i + 8, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum3= _mm256_add_pd(sum3, vp);

            vx = _mm256_load_pd(x + i + 12);
            vy = _mm256_load_pd(y + i + 12);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i + 12, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum4 = _mm256_add_pd(sum4, vp);

            vx = _mm256_load_pd(x + i + 16);
            vy = _mm256_load_pd(y + i + 16);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i + 16, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum5 = _mm256_add_pd(sum5, vp);

            vx = _mm256_load_pd(x + i + 20);
            vy = _mm256_load_pd(y + i + 20);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i + 20, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum6 = _mm256_add_pd(sum6, vp);

            vx = _mm256_load_pd(x + i + 24);
            vy = _mm256_load_pd(y + i + 24);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i + 24, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum7 = _mm256_add_pd(sum7, vp);

            vx = _mm256_loadu_pd(x + i + 28);
            vy = _mm256_loadu_pd(y + i + 28);
            vp = _mm256_fmadd_pd(va, vx, vy);
            _mm256_store_pd(res + i + 28, vp);
            vp = _mm256_andnot_pd(msk, vp); // fabs(vp)
            sum8 = _mm256_add_pd(sum8, vp);
        }

        sum1 = _mm256_add_pd(sum1, sum2);
        sum3 = _mm256_add_pd(sum3, sum4);
        sum5 = _mm256_add_pd(sum5, sum6);
        sum7 = _mm256_add_pd(sum7, sum8);

        sum1 = _mm256_add_pd(sum1, sum3);
        sum5 = _mm256_add_pd(sum5, sum7);

        sum1 = _mm256_add_pd(sum1, sum5);

        // Add up the horizontal sum
        s = sum1[0];
        for (i = 1; i < 4; i++) {
            s += sum1[i];
        }
    }

    // Handle the trailing entries
    for (i = n_32; i < n; i++) {
        res[i] = fma(a, x[i], y[i]);
        s += fabs(res[i]);
    }

    return s;
}
