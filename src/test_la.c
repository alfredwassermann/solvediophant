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

#include "const.h"
#include "arith.h"

DOUBLE *getArray(int len, int type) {
    DOUBLE *p;
    int i, sgn;

    // p = (DOUBLE *)calloc(len, sizeof(DOUBLE));
    p = (DOUBLE *)malloc(len * sizeof(DOUBLE));
    // p = (DOUBLE *)aligned_alloc(32, len * sizeof(DOUBLE));
    for (i = 0, sgn = 1.0; i < len; i++) {
        switch (type) {
            case 1:
                p[i] = 1.0 * (DOUBLE)(sgn * i);
                break;
            case 2:
                p[i] = sgn * 2.0 / (DOUBLE)(i + 1);
                break;
            case 3:
                p[i] = -1.0 / ((i+1) * (i+1));
                break;
            default:
                p[i] = 1.0;
        }
        sgn *= (-1.0);
        // printf("%lf ", p[i]);
    }
    return p;
}

//
// Naive implementations for comparison
//
DOUBLE sumNaive(DOUBLE *p, int n)
{
    DOUBLE s;
    int i;

    for (i = 0, s = 0.0; i < n; i++)
    {
        s += p[i];
    }
    return s;
}

/**
 * Naive summation using AVX
 */
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

DOUBLE dotNaive(DOUBLE *x, DOUBLE *y, int n)
{
    DOUBLE s;
    int i;
    if (n <= 0)
        return 0.0;

    for (i = 0, s = 0.0; i < n; i++)
    {
        s += x[i] * y[i];
    }
    return s;
}

DOUBLE dotNaiveQP(DOUBLE *x, DOUBLE *y, int n)
{
    _Float128 s;
    int i;
    if (n <= 0)
        return 0.0;

    for (i = 0, s = 0.0; i < n; i++)
    {
        s += x[i] * y[i];
    }
    return (DOUBLE)s;
}

int main(int argc, char *argv[])
{

    DOUBLE x = 0.00000000001;
    DOUBLE y = 100000.0;
    DOUBLE z;
    hiprec a;
    int n = 50000;


    // Test primitive functions I 
    #if 0
        printf("--------- ADD\n");
        z = x + y;
        printf("%0.20lf\n", z);
        a = twoSum(x, y);
        printf("%0.20lf %0.20lf\n", a.hi, a.lo);
        printf("--------- Split\n");
        a = split(z);
        printf("%0.20lf %0.20lf\n", a.hi, a.lo);
    #endif

    // Test primitive functions II
    #if 0
        printf("fma: %lf\n", fma(3.0, 1700000000.0, -1.0));
        printf("--------- TwoProduct\n");
        x = 11111111.111111111;
        y = 7.777777777;
        a = twoProd(x, y);
        printf("twoProd : %0.20lf %0.20lf\n", a.hi, a.lo);
        DOUBLE x1, y1;
        twoProd2(x, y, &x1, &y1);
        printf("twoProd2: %0.20lf %0.20lf\n", x1, y1);

        z = x * y;
        printf("FMA: %0.20lf %0.20lf\n", z, x * y - z);

        x = 10000000000000.0;
        printf("Naive:      %0.20lf\n", x * x);
        a = twoSquare(x);
        printf("twoSquare : %0.20lf %0.20lf\n", a.hi, a.lo);
        twoSquare2(x, &x1, &y1);
        printf("twoSquare2: %0.20lf %0.20lf\n", x1, y1);
        twoProd2(x, x, &x1, &y1);
        printf("twoSquare2: %0.20lf %0.20lf\n", x1, y1);

        printf("sqrt   : %0.20lf\n", sqrt(x1));        if (i >= n) {
            i -= 3;
            sum[0] += 
        }

        printf("accSqrt: %0.20lf\n", hiprec_sqrt(x1, y1));
    #endif

    // Test sqrt function
    #if 0
        printf("--------- sqrt\n");        if (i >= n) {
            i -= 3;
            sum[0] += 
        }

        y = 5000000000000000.0; 
        z = 2;
        printf("Hi1  : %0.20lf\n", hiprec_sqrt(y*y + 2* y * z, z));
        printf("Hi2  : %0.20lf\n", hiprec_sqrt(y*y, 2* y * z + z));
    #endif

    // Test summation
    #if 0
    {
        DOUBLE *p = getArray(n, 2);

        printf("--------- Sum2s\n");
        printf("Start\n");

        printf("Naive  %0.20lf\n", sumNaive(p, n));
        printf("NaiveA %0.20lf\n", sumNaiveAVX(p, n));
        printf("sum2s  %0.20lf\n", hiprec_sum2(p, n));
        printf("SUM    %0.20lf\n", hiprec_SUM(p, n));
        printf("SUMAVX %0.20lf\n", hiprec_sum_AVX(p, n));
        // printf("sum2  %0.20lf\n", hiprec_sumK(p, n, 2));
        // printf("sum3  %0.20lf\n", hiprec_sumK(p, n, 3));
        // printf("sum4  %0.20lf\n", hiprec_sumK(p, n, 4));
        // printf("sum5  %0.20lf\n", hiprec_sumK(p, n, 5));
    }
    #endif

    // Multiple summations
    #if 0
    {
        DOUBLE *p = getArray(50000, 2);

        printf("--------- Sum2s\n");
        printf("Start\n");
        double s = 0.0;

        for (int j = 0; j < 100000; j++) {
            // double q1 = sumNaive(p, n - j);
            // double q1 = sumNaiveAVX(p, n - j);
            double q1 = hiprec_sum_AVX(p, n - j);
            // double q1 = hiprec_sum2(p, n - j);
            // printf("%0.16lf %0.16lf  %0.16lf\n", q1, q2, q1 - q2);
            s += q1;
        }
        printf("Res: %0.20lf\n", s);
    }
    #endif

    // Test l1 norm
    #if 0
    {
        DOUBLE *p = getArray(n, 2);

        printf("hiprec l1  %0.20lf\n", hiprec_norm_l1(p, n));
        printf("hiprec l1  %0.20lf\n", hiprec_norm_l1_AVX(p, n));

        double s = 0.0;
        for (int j = 0; j < 100000; j++) {
            // double q1 = hiprec_norm_l1(p, n - j);
            double q1 = hiprec_norm_l1_AVX(p, n - j);
            s += q1;
        }
        printf("Res: %0.20lf\n", s);

    }
    #endif

    // Test dot product
    #if 1
    {
        printf("--------- Dot\n");
        const int n = 60013;

        DOUBLE q1;
        DOUBLE *p = getArray(n, 2);
        DOUBLE *q = getArray(n, 3);

        printf("Naive  %0.20lf\n", dotNaive(p, q, n));
        printf("NaiQP  %0.20lf\n", dotNaiveQP(p, q, n));
        printf("dot2   %0.20lf\n", hiprec_dot2(p, q, n));
        printf("dotAVX %0.20lf\n", hiprec_dot2_AVX(p, q, n));

        double s = 0.0;
        for (int j = 0; j < n; j++) {
            // q1 = hiprec_dot2(p, q, n - j);
            q1 = hiprec_dot2_AVX(p, q, n - j);
            s += q1;
        }
        printf("Multiple dot %0.20lf\n", s);
    }
    #endif

    // Test 2-norm
    #if 1
    {
        printf("--------- Norm\n");
        const int n = 60013;

        DOUBLE q1 = 0.0;
        DOUBLE *p = getArray(n, 2);

        printf("normsq      %0.20lf\n", hiprec_normsq_l2(p, n));
        printf("normsqAVX   %0.20lf\n", hiprec_normsq_l2_AVX(p, n));
        printf("\n");
        printf("NaiveNorm   %0.20lf\n", sqrt(dotNaive(p, p, n)));
        printf("NaiveNormQP %0.20lf\n", sqrt(dotNaiveQP(p, p, n)));
        printf("norm2       %0.20lf\n", hiprec_norm_l2(p, n));
        printf("normAVX     %0.20lf\n", hiprec_norm_l2_AVX(p, n));

        double s = 0.0;
        for (int j = 0; j < n; j++) {
            // q1 = hiprec_norm_l2(p, n - j);
            q1 = hiprec_norm_l2_AVX(p, n - j);
            s += q1;
        }
        printf("Multiple norm %0.20lf\n", s);
    }
    #endif

    return 0;
}
