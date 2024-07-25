#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <x86intrin.h>

#include "const.h"
#include "arith.h"

/**
* Ogita, Rump, Oishi: Accurate sum and dot product (2005)
*/

/**
 * Add two doubles and return the hiprec result.
 */
hiprec twoSum(DOUBLE a, DOUBLE b) {
    hiprec ret;
    DOUBLE z;

    ret.hi = a + b;
    z = ret.hi - a;
    ret.lo = (a - (ret.hi - z)) + (b - z);

    return ret;
}

/**
 * Add two doubles and return the hiprec result in parameters hi and lo.
 */
void twoSum2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo) {
    DOUBLE z;

    (*hi) = a + b;
    z = (*hi) - a;
    (*lo) = (a - ((*hi) - z)) + (b - z);
}

/**
 * Add two doubles a and b with |a| > |b| and return the hiprec result in parameters hi and lo.
 */
void fastTwoSum2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo) {
    // FastTwoSum for |a| > |b|
    (*hi) = a + b;
    (*lo) = (a - (*hi)) + b;
}

/**
 * Add two doubles a and b return hi in a and lo in b.
 */
void twoSum2i(DOUBLE *a, DOUBLE *b) {
    /** In-place summation */
    DOUBLE z, x;

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
hiprec split(DOUBLE a) {
    hiprec ret;
    const DOUBLE factor = 134217729;
    DOUBLE c;

    c = factor * a;
    ret.hi = (c - (c - a));
    ret.lo = (a - ret.hi);
    return ret;
}

/**
 * Multiply two numbers and return hiprec.
 */
hiprec twoProd(DOUBLE a, DOUBLE b) {
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
void twoProd2(DOUBLE a, DOUBLE b, DOUBLE *hi, DOUBLE *lo) {
    hiprec a_e, b_e;
    (*hi) = a * b;
    a_e = split(a);
    b_e = split(b);
    (*lo) = a_e.lo * b_e.lo - ((((*hi) - a_e.hi * b_e.hi) - a_e.lo * b_e.hi) - a_e.hi * b_e.lo);
}

/**
 * Multiply a * a and return hiprec result.
*/
hiprec twoSquare(DOUBLE a) {
    hiprec ret, a_e;
    ret.hi = a * a;
    a_e = split(a);
    ret.lo = a_e.lo * a_e.lo - ((ret.hi - a_e.hi * a_e.hi) - 2 * a_e.lo * a_e.hi);
    return ret;
}

/**
 * Multiply a * a and return hiprec result in parameters hi and lo.
 */
void twoSquare2(DOUBLE a, DOUBLE *hi, DOUBLE *lo) {
    hiprec a_e;
    (*hi) = a * a;
    a_e = split(a);
    (*lo) = a_e.lo * a_e.lo - (((*hi) - a_e.hi * a_e.hi) - 2 * a_e.lo * a_e.hi);
}

/**
 *  Accurate square root of hiprec number (T, t).
 */
DOUBLE hiprec_sqrt(DOUBLE T, DOUBLE t) {
    DOUBLE P, p, H, h;
    DOUBLE r;

    P = sqrt(T);

    twoSquare2(P, &H, &h);
    r = (T - H) - h;
    r = t + r;
    p = r / (2 * P);

    return P + p;
}

/**
 * Sum entries of array p of length n with double precision.
 */
DOUBLE hiprec_sum2(DOUBLE* p, int n) {
    DOUBLE sigma;
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
DOUBLE hiprec_SUM(DOUBLE* p, int n) {
    DOUBLE sigma;
    DOUBLE h, z;
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

//
// AVX2
//
DOUBLE sumNaiveAVX(DOUBLE *p, int n)
{
    long i = 0;
    long n_4 = n & -4;
    DOUBLE s;

    __m256d sum = {0.0, 0.0, 0.0, 0.0};
    while (i < n_4)
    {
        __m256d vp = _mm256_load_pd(&p[i]);
        sum = _mm256_add_pd(vp, sum);
        i += 4;
    }

    // sum = {a,b,c,d}
    __m128d xmm = _mm256_extractf128_pd(sum, 1); // xmm = {c,d}
    __m256d ymm = {xmm[0], xmm[1], 0, 0};        // ymm = {c,d,0,0}
    sum = _mm256_hadd_pd(sum, ymm);              // sum = {a+b,c+d,c+d,0}
    sum = _mm256_hadd_pd(sum, sum);              // sum = {a+b+c+d,a+b+c+d,c+d,c+d}

    s = sum[0];

    while (i < n)
    {
        s += p[i];
        i++;
    }

    return s;
}

DOUBLE hiprec_SUM_AVX(DOUBLE* p, int n) {
    hiprec s; 
    DOUBLE sigma_d;
    long i = 0;
    __m256d lo, h, z;

    if (n <= 0) return 0.0;

    long n_4 = n & -4;

    if (n_4 >= 4) {
        __m256d sum_hi = {0.0, 0.0, 0.0, 0.0};
        __m256d sigma  = {0.0, 0.0, 0.0, 0.0};

        for (i = 0; i < n_4; i += 4) {
            __m256d vp = _mm256_load_pd(&p[i]);
            TWOSUM_AVX(sum_hi, vp, sum_hi, lo, h, z);
            sigma = _mm256_add_pd(sigma, lo);
        }

        // Add the sums horizontally
        s.hi = sum_hi[0];
        sigma_d = sigma[0];
        for (i = 1; i < 4; i++) {
            s = twoSum(s.hi, sum_hi[i]);
            sigma_d += s.lo + sigma[i];
        }
    } else {
        s.hi = p[0];
        sigma_d = 0;
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
DOUBLE hiprec_sumK(DOUBLE* p, int n, int K) {
    DOUBLE s, alpha;
    DOUBLE q[K];
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
 * Sum absolute values of entries of array p of length n with double precision.
 */
DOUBLE hiprec_norm_l1(DOUBLE* p, int n) {
    DOUBLE sigma;
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
DOUBLE hiprec_normK_l1(DOUBLE* p, int n, int K) {
    DOUBLE s, alpha;
    DOUBLE q[K - 1];
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
 * Dot product of arrays x and y of length n with high precision.
*/
DOUBLE hiprec_dot2(DOUBLE* x, DOUBLE* y, int n) {
    DOUBLE p, q, s, h, r;
    int i;
    if (n <= 0) return 0.0;

    twoProd2(x[0], y[0], &p, &s);
    for (i = 1; i < n; i++) {
        twoProd2(x[i], y[i], &h, &r);
        twoSum2(p, h, &p, &q);
        s += (q + r);
    }
    return p + s;
}

/**
 * Dot product of arrays x and y of length n with K-fold precision.
*/
DOUBLE hiprec_dotK(DOUBLE* x, DOUBLE* y, int n, int K) {
    DOUBLE P, H;
    DOUBLE r[2 * n];
    int i;
    if (n <= 0) return 0.0;

    twoProd2(x[0], y[0], &P, &(r[0]));
    for (i = 1; i < n; i++) {
        twoProd2(x[i], y[i], &H, &(r[i]));
        twoSum2(P, H, &P, &(r[n + i - 1]));
    }
    r[2 * n - 1] = P;
    return hiprec_sumK(r, 2 * n, K);
}

/**
 * Dot product of arrays x and y of length n and step widths dx and dy with double precision.
 */
DOUBLE hiprec_dot2_row(DOUBLE* x, int dx, DOUBLE* y, int dy, int n) {
    DOUBLE p, q, s, h, r;
    int i, jx, jy;
    if (n <= 0) return 0.0;

    twoProd2(x[0], y[0], &p, &s);
    for (i = 1, jx = dx, jy = dy; i < n; i++, jx += dx, jy += dy) {
        twoProd2(x[jx], y[jy], &h, &r);
        twoSum2(p, h, &p, &q);
        s += (q + r);
    }
    return p + s;
}

/**
 * Square of ||x||_2, i.e.
 * dot product of array x with itself of length n with high precision.
 */
DOUBLE hiprec_normsq_l2(DOUBLE* x, int n) {
    DOUBLE S, s, P, p, H, h;
    DOUBLE c, d;
    int i;

    if (n <= 0) return 0.0;

    S = 0.0;
    s = 0.0;
    for (i = 0; i < n; i++) {
        twoSquare2(x[i], &P, &p);
        twoSum2(S, P, &H, &h);
        c = s + p;
        d = h + c;
        fastTwoSum2(H, d, &S, &s);
    }
    return S + s;
}

/**
 * ||x||_2, i.e. square root of dot product of array x with itself of length n with high precision.
 */
DOUBLE hiprec_norm_l2(DOUBLE* x, int n) {
    DOUBLE S, s, P, p, H, h;
    DOUBLE c, d;
    int i;

    if (n <= 0) return 0.0;

    S = 0.0;
    s = 0.0;
    for (i = 0; i < n; i++) {
        twoSquare2(x[i], &P, &p);
        twoSum2(S, P, &H, &h);
        c = s + p;
        d = h + c;
        fastTwoSum2(H, d, &S, &s);
    }
    return hiprec_sqrt(S, s);
}

/**
 * ||x||_2, i.e. square root of dot product of array x with itself of length n with K-fold precision.
 */
DOUBLE hiprec_normK_l2(DOUBLE* x, int n, int K) {
    DOUBLE S, P;
    DOUBLE r[2 * n];
    int i;

    if (n <= 0) return 0.0;

    S = 0.0;
    for (i = 0; i < n; i++) {
        twoSquare2(x[i], &P, &(r[i]));
        twoSum2(S, P, &S, &(r[n + i - 1]));
    }
    r[2 * n - 1] = S;

    return hiprec_sqrt(hiprec_sumK(r, 2 * n, K), 0.0);
}

