#include <stdio.h>
#include <stdlib.h>

#include "const.h"

typedef struct dblexact {
    DOUBLE x;
    DOUBLE y;
} doubleExact;


/*
* Ogita, Rump, Oishi:
* Accurate sum and dot product (205)
*/

doubleExact twoSum(DOUBLE a, DOUBLE b) {
    doubleExact ret;
    DOUBLE z;

    ret.x = a + b;

    z = ret.x - a;

    ret.y = (a - (ret.x - z)) + (b - z);
    // FastTwoSum for a > b
    // ret.y = (a - ret.x) + b;

    return ret;
}

void twoSum2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y) {
    DOUBLE z;

    (*x) = a + b;

    z = (*x) - a;

    (*y) = (a - ((*x) - z)) + (b - z);
    // FastTwoSum for a > b
    // (*y) = (a - (*x)) + b;
}

/** in place summation*/
void twoSum2i(DOUBLE *a, DOUBLE *b) {
    DOUBLE z, x;

    x = (*a) + (*b);

    z = x - (*a);
    (*b) = ((*a) - (x - z)) + ((*b) - z);
    // FastTwoSum for a > b
    // (*b) = ((*a) - x) + (*b);

    (*a) = x;
}

doubleExact split(DOUBLE a) {
    doubleExact ret;
    const DOUBLE factor = 134217729;
    DOUBLE c;

    c = factor * a;
    ret.x = (c - (c - a));
    ret.y = (a - ret.x);
    return ret;
}

doubleExact twoProd(DOUBLE a, DOUBLE b) {
    doubleExact ret, a_e, b_e;
    ret.x = a * b;
    a_e = split(a);
    b_e = split(b);
    ret.y = a_e.y * b_e.y - (((ret.x - a_e.x * b_e.x) - a_e.y * b_e.x) - a_e.x * b_e.y);
    return ret;
}

void twoProd2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y) {
    doubleExact a_e, b_e;
    (*x) = a * b;
    a_e = split(a);
    b_e = split(b);
    (*y) = a_e.y * b_e.y - ((((*x) - a_e.x * b_e.x) - a_e.y * b_e.x) - a_e.x * b_e.y);
}

DOUBLE sumNaive(DOUBLE* p, int n) {
    DOUBLE s;
    int i;

    for (i = 0, s = 0.0; i < n; i++) {
        s += p[i];
    }
    return s;
}

DOUBLE sum2s(DOUBLE* p, int n) {
    DOUBLE sigma;
    doubleExact s;
    int i;

    if (n <= 0) return 0.0;

    s.x = p[0];
    sigma = 0.0;

    for (i = 1; i < n; i++) {
        s = twoSum(s.x, p[i]);
        sigma += s.y;
    }

    return s.x + sigma;
}

DOUBLE sumKvert(DOUBLE* p, int n, int K) {
    DOUBLE s, alpha;
    DOUBLE q[K - 1];
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

DOUBLE dot2(DOUBLE* x, DOUBLE* y, int n) {
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

int main(int argc, char *argv[]) {

    DOUBLE x = 0.00000000001;
    DOUBLE y = 100000.0;
    DOUBLE z, sgn;
    doubleExact a;

    #if 0
        printf("--------- ADD\n");
        z = x + y;
        printf("%0.20lf\n", z);
        a = twoSum(x, y);
        printf("%0.20lf %0.20lf\n", a.x, a.y);

        printf("--------- Split\n");
        a = split(z);
        printf("%0.20lf %0.20lf\n", a.x, a.y);
    #endif

    #if 0
        printf("--------- TwoProduct\n");
        x = 11111111.111111111;
        y = 7.777777777;
        a = twoProd(x, y);
        printf("twoProd : %0.20lf %0.20lf\n", a.x, a.y);
        DOUBLE x1, y1;
        twoProd2(x, y, &x1, &y1);
        printf("twoProd2: %0.20lf %0.20lf\n", x1, y1);

        z = x * y;
        printf("FMA: %0.20lf %0.20lf\n", z, x * y - z);
    #endif

    #if 0
        printf("--------- Sum2s\n");
        int i;
        const int n = 10000;
        DOUBLE p[n];

        for (i = 0, sgn = 1.0; i < n; i++) {
            // p[i] = 2.0 / (DOUBLE)(i + 1);
            // p[i] = 10000.0 * (DOUBLE)(sgn * i);
            p[i] = 1.0 / ((i+1) * (i+1));
            // sgn *= (-1.0);
            // printf("%lf ", p[i]);
        }
        // printf("\n");

        printf("Naive %0.20lf\n", sumNaive(p, n));
        printf("sum2s %0.20lf\n", sum2s(p, n));
        printf("sum2  %0.20lf\n", sumKvert(p, n, 2));
        printf("sum3  %0.20lf\n", sumKvert(p, n, 3));
        printf("sum4  %0.20lf\n", sumKvert(p, n, 4));
        printf("sum5  %0.20lf\n", sumKvert(p, n, 5));
    #endif

    #if 1
        printf("--------- Dot\n");
        int i;
        const int n = 60000;
        DOUBLE p[n];
        DOUBLE q[n];

        for (i = 0, sgn = 1.0; i < n; i++) {
            p[i] = 2.0 / (DOUBLE)(i + 1);
            q[i] = 1.0 / ((i+1) * (i+1));
        }

        printf("Naive %0.20lf\n", dotNaive(p, q, n));
        printf("NaiQP %0.20lf\n", dotNaiveQP(p, q, n));
        printf("dot2  %0.20lf\n", dot2(p, q, n));
    #endif

    return 0;
}