// gcc -O3 -march=skylake-avx512 -DTEST_MAIN sum_avx512_u4_skx.c -o sum_u4_skx

// sum_avx512_u4_skx.c
#include <immintrin.h>
#include <stddef.h>

static inline double hsum512_pd(__m512d v) {
    __m256d lo256 = _mm512_castpd512_pd256(v);
    __m256d hi256 = _mm512_extractf64x4_pd(v, 1);
    __m256d s256 = _mm256_add_pd(lo256, hi256);
    __m128d lo128 = _mm256_castpd256_pd128(s256);
    __m128d hi128 = _mm256_extractf128_pd(s256, 1);
    __m128d v2 = _mm_add_pd(lo128, hi128);
    __m128d sh = _mm_shuffle_pd(v2, v2, 0x1);
    return _mm_cvtsd_f64(_mm_add_sd(v2, sh));
}

double sum_avx512_u4_skx(const double *a, size_t n) {
#if defined(__AVX512F__)
    __m512d s0 = _mm512_setzero_pd();
    __m512d s1 = _mm512_setzero_pd();
    __m512d s2 = _mm512_setzero_pd();
    __m512d s3 = _mm512_setzero_pd();

    size_t i = 0;
    size_t n32 = n & ~(size_t)31;

    // 4× unrolled main loop: 32 doubles per iter
    for (; i < n32; i += 32) {
        s0 = _mm512_add_pd(s0, _mm512_loadu_pd(a + i +  0));
        s1 = _mm512_add_pd(s1, _mm512_loadu_pd(a + i +  8));
        s2 = _mm512_add_pd(s2, _mm512_loadu_pd(a + i + 16));
        s3 = _mm512_add_pd(s3, _mm512_loadu_pd(a + i + 24));
    }

    // Handle leftover full vectors (chunks of 8)
    size_t n8 = n & ~(size_t)7;
    for (; i < n8; i += 8) {
        s0 = _mm512_add_pd(s0, _mm512_loadu_pd(a + i));
    }

    // Final masked tail (<8)
    unsigned r = (unsigned)(n - n8);
    if (r) {
        __mmask8 k = (1u << r) - 1u;
        s0 = _mm512_add_pd(s0, _mm512_maskz_loadu_pd(k, a + n8));
    }

    __m512d s01 = _mm512_add_pd(s0, s1);
    __m512d s23 = _mm512_add_pd(s2, s3);
    return hsum512_pd(_mm512_add_pd(s01, s23));
#else
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += a[i];
    return s;
#endif
}

#ifdef TEST_MAIN
#include <stdio.h>
int main(void) {
    enum { N = 100003 };
    static double a[N];
    for (int i = 0; i < N; ++i) a[i] = (double)(i & 1023) * 0.25;
    printf("%.17g\n", sum_avx512_u4_skx(a, N));
    return 0;
}
#endif