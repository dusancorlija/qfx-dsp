/*
 * noise_shaping.c - measures the output error of a Q31 biquad against a
 * double-precision reference, with and without error feedback, and prints
 * the error spectrum as CSV (freq_hz, round_db, shaped_db).
 *
 *   ./noise_shaping > docs/noise_shaping.csv
 *   python3 tools/plot.py docs/noise_shaping.csv docs/noise_shaping.svg
 *
 * A summary table of SNR versus cutoff frequency goes to stderr.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qfx/biquad.h"
#include "qfx/design.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FS 48000.0
#define N 65536
#define SETTLE 8192

static q31_t x[N], y[N];
static double xr[N], yr[N];

static void fft(double *re, double *im, int n)
{
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        const double ang = -2 * M_PI / len;
        for (int i = 0; i < n; i += len)
            for (int k = 0; k < len / 2; k++) {
                const double wr = cos(ang * k), wi = sin(ang * k);
                const double ur = re[i + k], ui = im[i + k];
                const double vr = re[i + k + len / 2] * wr - im[i + k + len / 2] * wi;
                const double vi = re[i + k + len / 2] * wi + im[i + k + len / 2] * wr;
                re[i + k] = ur + vr; im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr; im[i + k + len / 2] = ui - vi;
            }
    }
}

/* Runs one configuration; returns SNR in dB and optionally the error.
 * With exact_coef the reference uses the ideal double coefficients, so the
 * result also includes the error from quantising the coefficients to Q31.
 * Without it the reference uses the quantised coefficients and the result
 * is pure arithmetic (rounding) noise. */
static double run(qfx_sos_t ideal, uint8_t mode, int exact_coef, double *err)
{
    qfx_biquad_coef_t c;
    qfx_biquad_state_t st;
    qfx_biquad_t f;
    qfx_biquad_coef_from_double(&c, ideal.b0, ideal.b1, ideal.b2, ideal.a1, ideal.a2);
    qfx_biquad_init(&f, &c, &st, 1, mode);
    qfx_biquad_process(&f, x, y, N);

    const double k = ldexp(1.0, -(31 - c.shift));
    const qfx_sos_t quant = { c.b0 * k, c.b1 * k, c.b2 * k, c.a1 * k, c.a2 * k };
    const qfx_sos_t s = exact_coef ? ideal : quant;

    double x1 = 0, x2 = 0, y1 = 0, y2 = 0, ps = 0, pe = 0;
    for (int i = 0; i < N; i++) {
        const double x0 = xr[i];
        const double y0 = s.b0 * x0 + s.b1 * x1 + s.b2 * x2 - s.a1 * y1 - s.a2 * y2;
        x2 = x1; x1 = x0; y2 = y1; y1 = y0;
        yr[i] = y0;
        const double e = q31_to_double(y[i]) - y0;
        if (err) err[i] = e;
        if (i >= SETTLE) { ps += y0 * y0; pe += e * e; }
    }
    return 10 * log10(ps / pe);
}

static void make_input(double f0)
{
    srand(3);
    for (int i = 0; i < N; i++) {
        /* tone + small dither so the error is noise-like, not periodic */
        const double d = ((double)rand() / RAND_MAX - 0.5) * 1e-6;
        const double v = 0.5 * sin(2 * M_PI * f0 * i / FS) + d;
        x[i] = q31_from_double(v);
        xr[i] = q31_to_double(x[i]);
    }
}

int main(void)
{
    /* SNR table vs cutoff */
    const double cut[] = { 10, 30, 100, 300, 1000, 3000, 10000 };
    fprintf(stderr, "| Low-pass cutoff (fs = 48 kHz) | Rounded | Error feedback | Gain | Incl. coefficient quantisation |\n");
    fprintf(stderr, "|---:|---:|---:|---:|---:|\n");
    for (size_t k = 0; k < sizeof(cut) / sizeof(cut[0]); k++) {
        qfx_sos_t s = qfx_design_lowpass(FS, cut[k], M_SQRT1_2);
        make_input(cut[k] / 4);
        const double a = run(s, QFX_BIQUAD_ROUND, 0, NULL);
        const double b = run(s, QFX_BIQUAD_NOISE_SHAPING, 0, NULL);
        const double q = run(s, QFX_BIQUAD_NOISE_SHAPING, 1, NULL);
        fprintf(stderr, "| %.0f Hz | %.1f dB | %.1f dB | %+.1f dB | %.1f dB |\n", cut[k], a, b, b - a, q);
    }

    /* Error spectrum for a 30 Hz low-pass */
    static double e0[N], e1[N], re[N], im[N];
    qfx_sos_t s = qfx_design_lowpass(FS, 30, M_SQRT1_2);
    make_input(7.5);
    run(s, QFX_BIQUAD_ROUND, 0, e0);
    run(s, QFX_BIQUAD_NOISE_SHAPING, 0, e1);

    enum { M = 16384 };
    double psd[2][M / 2] = { { 0 } };
    int segs = 0;
    for (int off = SETTLE; off + M <= N; off += M / 2, segs++)
        for (int v = 0; v < 2; v++) {
            const double *e = v ? e1 : e0;
            for (int i = 0; i < M; i++) {
                const double w = 0.5 - 0.5 * cos(2 * M_PI * i / (M - 1));
                re[i] = e[off + i] * w;
                im[i] = 0;
            }
            fft(re, im, M);
            for (int i = 0; i < M / 2; i++) psd[v][i] += re[i] * re[i] + im[i] * im[i];
        }

    printf("freq_hz,round_db,shaped_db\n");
    /* log-spaced bins keep the CSV small */
    double last = 0;
    for (int i = 1; i < M / 2; i++) {
        const double f = i * FS / M;
        if (f < last * 1.03) continue;
        last = f;
        const double a = 10 * log10(psd[0][i] / segs + 1e-300);
        const double b = 10 * log10(psd[1][i] / segs + 1e-300);
        printf("%.2f,%.2f,%.2f\n", f, a, b);
    }
    return 0;
}
