/*
 * bench.c - host throughput comparison, float vs fixed point.
 *
 * Host numbers only show relative cost of the algorithms; on a Cortex-M
 * without an FPU (M0+/M3) the float rows become 10-50x slower, while on an
 * M4F/M7 float and Q31 are close and Q15 FIR wins through SIMD MACs.
 */
#define _POSIX_C_SOURCE 199309L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "qfx/biquad.h"
#include "qfx/design.h"
#include "qfx/fir.h"

#define FS 48000.0
#define N 4096
#define REPS 400

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static volatile float sinkf;
static volatile q31_t sink31;
static volatile q15_t sink15;

static void float_cascade(const qfx_sos_t *s, int stages, float *st, float *x, int n)
{
    for (int k = 0; k < stages; k++) {
        const float b0 = (float)s[k].b0, b1 = (float)s[k].b1, b2 = (float)s[k].b2;
        const float a1 = (float)s[k].a1, a2 = (float)s[k].a2;
        float x1 = st[4 * k], x2 = st[4 * k + 1], y1 = st[4 * k + 2], y2 = st[4 * k + 3];
        for (int i = 0; i < n; i++) {
            const float x0 = x[i];
            const float y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
            x2 = x1; x1 = x0; y2 = y1; y1 = y0;
            x[i] = y0;
        }
        st[4 * k] = x1; st[4 * k + 1] = x2; st[4 * k + 2] = y1; st[4 * k + 3] = y2;
    }
}

static void float_fir(const float *h, float *st, int taps, int *idx, const float *x, float *y, int n)
{
    for (int i = 0; i < n; i++) {
        *idx = (*idx == 0) ? taps - 1 : *idx - 1;
        st[*idx] = st[*idx + taps] = x[i];
        float acc = 0;
        for (int k = 0; k < taps; k++) acc += h[k] * st[*idx + k];
        y[i] = acc;
    }
}

int main(void)
{
    static float xf[N], yf[N];
    static q31_t x31[N];
    static q15_t x15[N], y15[N];
    srand(42);
    for (int i = 0; i < N; i++) {
        const double v = 0.5 * sin(0.01 * i) + 0.1 * ((double)rand() / RAND_MAX - 0.5);
        xf[i] = (float)v;
        x31[i] = q31_from_double(v);
        x15[i] = q15_from_double(v);
    }

    /* 8th-order Butterworth low-pass = 4 biquads */
    qfx_sos_t sos[4];
    qfx_design_butter_lowpass(FS, 4000, 8, sos);
    qfx_biquad_coef_t c[4];
    qfx_biquad_state_t st[4];
    for (int k = 0; k < 4; k++)
        qfx_biquad_coef_from_double(&c[k], sos[k].b0, sos[k].b1, sos[k].b2, sos[k].a1, sos[k].a2);

    printf("| Kernel | ns / sample |\n|---|---:|\n");

    float fst[16] = { 0 };
    double t0 = now_s();
    for (int r = 0; r < REPS; r++) float_cascade(sos, 4, fst, xf, N);
    double t = now_s() - t0;
    sinkf = xf[N - 1];
    printf("| biquad x4, float | %.2f |\n", t / (REPS * (double)N) * 1e9);

    for (int m = 0; m < 2; m++) {
        qfx_biquad_t f;
        qfx_biquad_init(&f, c, st, 4, (uint8_t)m);
        t0 = now_s();
        for (int r = 0; r < REPS; r++) qfx_biquad_process(&f, x31, x31, N);
        t = now_s() - t0;
        sink31 = x31[N - 1];
        printf("| biquad x4, Q31 %s | %.2f |\n", m ? "error feedback" : "rounded",
               t / (REPS * (double)N) * 1e9);
    }

    enum { TAPS = 64 };
    float hf[TAPS], stf[2 * TAPS] = { 0 };
    q15_t h15[TAPS], st15[QFX_FIR_STATE_LEN(TAPS)];
    for (int k = 0; k < TAPS; k++) {
        const double w = 0.54 - 0.46 * cos(2 * 3.141592653589793 * k / (TAPS - 1));
        hf[k] = (float)(w / TAPS);
        h15[k] = q15_from_double(w / TAPS);
    }
    int idx = 0;
    t0 = now_s();
    for (int r = 0; r < REPS; r++) float_fir(hf, stf, TAPS, &idx, xf, yf, N);
    t = now_s() - t0;
    sinkf = yf[N - 1];
    printf("| FIR 64 taps, float | %.2f |\n", t / (REPS * (double)N) * 1e9);

    qfx_fir_t fir;
    qfx_fir_init(&fir, h15, st15, TAPS);
    t0 = now_s();
    for (int r = 0; r < REPS; r++) qfx_fir_process(&fir, x15, y15, N);
    t = now_s() - t0;
    sink15 = y15[N - 1];
    printf("| FIR 64 taps, Q15 | %.2f |\n", t / (REPS * (double)N) * 1e9);
    return 0;
}
