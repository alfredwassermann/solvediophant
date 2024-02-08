#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "const.h"
#include "linalg.h"

int main(int argc, char *argv[]) {

    DOUBLE x = 0.00000000001;
    DOUBLE y = 100000.0;
    DOUBLE z, sgn;
    doubleExact a;
    int i;

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

    #if 1
        printf("fma: %lf\n", fma(3.0, 1700000000.0, -1.0));
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

        x = 10000000000000.0;
        printf("Naive:      %0.20lf\n", x * x);
        a = twoSquare(x);
        printf("twoSquare : %0.20lf %0.20lf\n", a.x, a.y);
        twoSquare2(x, &x1, &y1);
        printf("twoSquare2: %0.20lf %0.20lf\n", x1, y1);
        twoProd2(x, x, &x1, &y1);
        printf("twoSquare2: %0.20lf %0.20lf\n", x1, y1);

        printf("sqrt   : %0.20lf\n", sqrt(x1));
        printf("accSqrt: %0.20lf\n", hiprec_sqrt(x1, y1));
    #endif

    #if 0
        printf("--------- sqrt\n");
        x = 5000000000000000.0 + 2.0;
        printf("Naive: %0.20lf\n", sqrt(x * x));
        y = 5000000000000000.0; 
        z = 2;
        printf("Hi1  : %0.20lf\n", hiprec_sqrt(y*y + 2* y * z, z));
        printf("Hi2  : %0.20lf\n", hiprec_sqrt(y*y, 2* y * z + z));
    #endif

    #if 0
    {
        printf("--------- Sum2s\n");
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
