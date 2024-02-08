#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "const.h"
#include "linalg.h"

/*
* Ogita, Rump, Oishi:
* Accurate sum and dot product (2005)
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
}

void fastTwoSum2(DOUBLE a, DOUBLE b, DOUBLE *x, DOUBLE *y) {
    // FastTwoSum for |a| > |b|
    (*x) = a + b;
    (*y) = (a - (*x)) + b;
}

/** In-place summation */
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

doubleExact twoSquare(DOUBLE a) {
    doubleExact ret, a_e;
    ret.x = a * a;
    a_e = split(a);
    ret.y = a_e.y * a_e.y - ((ret.x - a_e.x * a_e.x) - 2 * a_e.y * a_e.x);
    return ret;
}

void twoSquare2(DOUBLE a, DOUBLE *x, DOUBLE *y) {
    doubleExact a_e;
    (*x) = a * a;
    a_e = split(a);
    (*y) = a_e.y * a_e.y - (((*x) - a_e.x * a_e.x) - 2 * a_e.y * a_e.x);
}

DOUBLE hiprec_sqrt(DOUBLE T, DOUBLE t) {
    DOUBLE P, p, H, h;
    DOUBLE r;

    P = sqrt(T);

    twoSquare2(P, &H, &h);
    r = (T - H) - h;
    // printf("%0.20lf %0.20f %0.20lf %0.20lf\n", H, h, r, t);
    r = t + r;
    p = r / (2 * P);

    // printf("%0.20lf %0.20f %0.20lf\n", P, p, r);
    return P + p;
}

DOUBLE hiprec_sum2(DOUBLE* p, int n) {
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

DOUBLE hiprec_sumK(DOUBLE* p, int n, int K) {
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

DOUBLE hiprec_norm_l1(DOUBLE* p, int n) {
    DOUBLE sigma;
    doubleExact s;
    int i;

    if (n <= 0) return 0.0;

    s.x = abs(p[0]);
    sigma = 0.0;

    for (i = 1; i < n; i++) {
        s = twoSum(s.x, abs(p[i]));
        sigma += s.y;
    }

    return s.x + sigma;
}

DOUBLE hiprec_normK_l1(DOUBLE* p, int n, int K) {
    DOUBLE s, alpha;
    DOUBLE q[K - 1];
    int i, j, k;

    if (n <= 0) return 0.0;

    K = (n < K) ? n : K;

    for (i = 0; i < K; i++) {
        q[i] = 0.0;
    }

    for (i = 0; i < K - 1; i++) {
        s = abs(p[i]);
        for (k = 0; k < i - 1; k++) {
            twoSum2i(&q[k], &s);
        }
        q[i] = s;
    }
    s = q[K - 1];

    for (i = K - 1; i < n; i++) {
        alpha = abs(p[i]);
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

DOUBLE hiprec_dotK(DOUBLE* x, DOUBLE* y, int n, int K) {
    DOUBLE p, q, s, h;
    DOUBLE r[2 * n];
    int i;
    if (n <= 0) return 0.0;

    twoProd2(x[0], y[0], &p, &(r[0]));
    for (i = 1; i < n; i++) {
        twoProd2(x[i], y[i], &h, &(r[i]));
        twoSum2(p, h, &p, &(r[n + i - 1]));
    }
    r[2 * n - 1] = p;
    return hiprec_sumK(r, 2 * n, K);
}

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

DOUBLE hiprec_normK_l2(DOUBLE* x, int n, int K) {
    DOUBLE S, s, P, p, H, h;
    DOUBLE c, d;
    DOUBLE r[2 * n];
    int i;

    if (n <= 0) return 0.0;

    S = 0.0;
    s = 0.0;
    for (i = 0; i < n; i++) {
        twoSquare2(x[i], &P, &(r[i]));
        twoSum2(S, P, &S, &(r[n + i - 1]));
    }
    r[2 * n - 1] = S;
    
    return hiprec_sqrt(hiprec_sumK(r, 2 * n, K), 0.0);
}

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
