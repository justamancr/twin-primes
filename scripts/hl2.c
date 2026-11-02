// Exact second-order Hardy–Littlewood prediction of chi2/dof for twin-prime centers mod q, vs measurement.
// usage: ./hl2 x   (x up to ~1e9 on 7 GB)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

static const double C2 = 0.6601618158468696;
static const int QS[] = {31, 61, 113, 251, 503, 1009};
#define NQ 6

int main(int argc, char **argv) {
    uint64_t x = (uint64_t) atof(argv[1]);
    uint64_t nodd = x / 2 + 2;
    uint8_t *s = calloc(nodd / 8 + 1, 1);      // bit i -> 2i+1 composite?  (0 = prime candidate)
    #define ISP(n) (!(s[((n)>>1)>>3] & (1u << (((n)>>1)&7))))
    #define SETC(n) (s[((n)>>1)>>3] |= (1u << (((n)>>1)&7)))
    SETC(1);
    for (uint64_t p = 3; p * p <= x; p += 2)
        if (ISP(p)) for (uint64_t m = p * p; m <= x; m += 2 * p) SETC(m);
    // twins (lower member > 1000), histograms of centers mod q
    uint64_t T = 0; uint64_t hist[NQ][1024]; memset(hist, 0, sizeof hist);
    for (uint64_t p = 1001; p + 2 <= x; p += 2) if (ISP(p) && ISP(p + 2)) {
        T++; uint64_t c = p + 1;
        for (int i = 0; i < NQ; i++) hist[i][c % QS[i]]++;
    }
    fprintf(stderr, "x=%llu sieve+twins done, T=%llu\n", (unsigned long long)x, (unsigned long long)T);
    // rho(h), h = 6k, k = 0..K
    uint64_t K = x / 6;
    double *rho = malloc((K + 1) * sizeof(double));
    for (uint64_t k = 0; k <= K; k++) rho[k] = 6.0;
    long double glog = 0.0L;
    for (uint64_t p = 5; p <= x; p += 2) {
        if (!ISP(p)) continue;
        double pd = (double)p;
        double g = (1 - 4 / pd) / ((1 - 2 / pd) * (1 - 2 / pd));
        double m2 = (1 - 3 / pd) / ((1 - 2 / pd) * (1 - 2 / pd));
        double z0 = 1 / (1 - 2 / pd);
        glog += logl((long double)g);
        uint64_t inv6 = (p % 6 == 1) ? (5 * p + 1) / 6 : (p + 1) / 6;
        uint64_t r1 = (2 * inv6) % p, r2 = p - r1;
        double fm = m2 / g, fz = z0 / g;
        if (p <= K) {
            for (uint64_t k = p; k <= K; k += p) rho[k] *= fz;
            for (uint64_t k = r1; k <= K; k += p) rho[k] *= fm;
            for (uint64_t k = r2; k <= K; k += p) rho[k] *= fm;
        } else {
            if (r1 <= K) rho[r1] *= fm;
            if (r2 <= K) rho[r2] *= fm;
        }
    }
    double G = (double) expl(glog);
    long double mean = 0;
    for (uint64_t k = 1; k <= K; k++) { rho[k] *= G; mean += rho[k]; }
    rho[0] = 0;
    fprintf(stderr, "rho built; mean rho = %.7Lf\n", mean / (long double)x);
    // sums: Aall = 2 sum_k rho(6k)(x-6k);  Aq = 2 sum_{q|k} rho(6k)(x-6k)
    long double Aall = 0;
    for (uint64_t k = 1; k <= K; k++) Aall += (long double)rho[k] * (long double)(x - 6 * k);
    Aall *= 2;
    printf("x=%.0e  T=%llu\n", (double)x, (unsigned long long)T);
    printf("%6s %12s %10s\n", "q", "prediction", "measured");
    for (int i = 0; i < NQ; i++) {
        int q = QS[i]; long double kap = q - 2;
        long double Aq = 0;
        for (uint64_t k = q; k <= K; k += q) Aq += (long double)rho[k] * (long double)(x - 6 * k);
        Aq *= 2;
        long double Tx = (long double)T / (long double)x;
        long double EV = (long double)T * (1 - 1 / kap) + Tx * Tx * (Aq - Aall / kap);
        long double pred = EV / ((kap - 1) * (long double)T / kap);
        // measured
        long double e = 0; int n = 0;
        for (int c = 0; c < q; c++) if (c != 1 && c != q - 1) { e += hist[i][c]; n++; }
        e /= n; long double chi = 0;
        for (int c = 0; c < q; c++) if (c != 1 && c != q - 1) { long double d = hist[i][c] - e; chi += d * d / e; }
        printf("%6d %12.3Lf %10.3Lf\n", q, pred, chi / (n - 1));
    }
    return 0;
}
