#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <x86intrin.h>

#include "const.h"
#include "arith.h"

//
// Naive implementations for comparison
//
DOUBLE sumNaive(DOUBLE* p, int n) {
    DOUBLE s;
    int i;

    for (i = 0, s = 0.0; i < n; i++) {
        s += p[i];
    }
    return s;
}

DOUBLE dotNaive(DOUBLE* x, DOUBLE* y, int n) {
    DOUBLE s;
    int i;
    if (n <= 0) return 0.0;

    for (i = 0, s = 0.0; i < n; i++) {
        s += x[i] * y[i];
    }
    return s;
}

DOUBLE dotNaiveQP(DOUBLE* x, DOUBLE* y, int n) {
    _Float128 s;
    int i;
    if (n <= 0) return 0.0;

    for (i = 0, s = 0.0; i < n; i++) {
        s += x[i] * y[i];
    }
    return (DOUBLE)s;
}

//
// AVX2
//
DOUBLE sumNaiveAVX(DOUBLE* p, int n) {
       int i;
       DOUBLE s;
        
        __m256d sum = {0.0, 0.0, 0.0, 0.0};
        for (i = 0; i < n; i += 4) {
            __m256d vp = _mm256_load_pd(&p[i]);
            sum = _mm256_add_pd (vp, sum);
        }
        
        __m128d xmm = _mm256_extractf128_pd (sum, 1); 
        __m256d ymm = {xmm[0], xmm[1], 0, 0};
        sum = _mm256_hadd_pd (sum, ymm); 
        sum = _mm256_hadd_pd (sum, sum);

        s = sum[0];
        return s;
}


int main(int argc, char *argv[]) {

    DOUBLE x = 0.00000000001;
    DOUBLE y = 100000.0;
    DOUBLE z;
    hiprec a;
    int i;

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

    #if 1
    {
        printf("--------- Sum2s\n");
        const int n = 500000;
        // DOUBLE p[n];
        DOUBLE *p;
        int sgn;

        p = (DOUBLE*)calloc(n, sizeof(DOUBLE));
        for (i = 0, sgn = 1.0; i < n; i++) {
            p[i] = sgn * 2.0 / (DOUBLE)(i + 1);
            // p[i] = 1.0;
            // if (i % 31 == 0) {
            //     p[i] *= 1000000;
            // }
            // if (i % 17 == 0) {
            //     p[i] *= 25252525 * sgn;
            // }
            // p[i] = 1.0 * (DOUBLE)(sgn * i);
            // p[i] = 1.0 / ((i+1) * (i+1));
            sgn *= (-1.0);
            // printf("%lf ", p[i]);
        }
        // printf("\n");

        printf("Start\n");
        double s = 0.0;
        for (int j = 0; j < 10000; j++) {
            // double q1 = sumNaive(p, n - j);
            double q1 = sumNaiveAVX(p, n - j);
            // printf("%0.16lf %0.16lf  %0.16lf\n", q1, q2, q1 - q2);
            s += q1;
        }
        printf("Naive %0.20lf\n", sumNaive(p, n));
        printf("Naive %0.20lf\n", s);
        // printf("sum2s %0.20lf\n", hiprec_sum2(p, n));
        // printf("sum2  %0.20lf\n", hiprec_sumK(p, n, 2));
        // printf("sum3  %0.20lf\n", hiprec_sumK(p, n, 3));
        // printf("sum4  %0.20lf\n", hiprec_sumK(p, n, 4));
        // printf("sum5  %0.20lf\n", hiprec_sumK(p, n, 5));
    }
    #endif

    #if 0
    {
        printf("--------- Dot\n");
        const int n = 60000;
        DOUBLE p[n];
        DOUBLE q[n];

        for (i = 0, sgn = 1.0; i < n; i++) {
            p[i] = 2.0 / (DOUBLE)(i + 1);
            q[i] = 1.0 / ((i+1) * (i+1));
        }

        printf("Naive %0.20lf\n", dotNaive(p, q, n));
        printf("NaiQP %0.20lf\n", dotNaiveQP(p, q, n));
        printf("dot2  %0.20lf\n", hiprec_dot2(p, q, n));
        printf("NaiveNorm   %0.20lf\n", sqrt(dotNaive(p, p, n)));
        printf("NaiveNormQP %0.20lf\n", sqrt(dotNaiveQP(p, p, n)));
        printf("norm        %0.20lf\n", hiprec_norm_l2(p, n));
    }
    #endif

    return 0;
}
