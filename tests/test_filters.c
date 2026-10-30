#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "qfx/biquad.h"
#include "qfx/cic.h"
#include "qfx/design.h"
#include "qfx/fir.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FS 48000.0
#define N 8192

static q31_t in31[N], out31[N], tmp31[N];
static double ref[N];

/* Coefficients exactly as the Q31 filter sees them. Comparing against these
 * isolates arithmetic noise from coefficient quantisation. */
static qfx_sos_t quantised(qfx_sos_t s)
{
    qfx_biquad_coef_t c;
    qfx_biquad_coef_from_double(&c, s.b0, s.b1, s.b2, s.a1, s.a2);
    const double k = ldexp(1.0, -(31 - c.shift));
    qfx_sos_t q = { c.b0 * k, c.b1 * k, c.b2 * k, c.a1 * k, c.a2 * k };
    return q;
}

/* Double-precision DF1 reference using the quantised coefficients. */
static void ref_cascade(const qfx_sos_t *sos, int stages, const double *x, double *y, int n)
{
    qfx_sos_t s[4];
    for (int k = 0; k < stages; k++) s[k] = quantised(sos[k]);
    memcpy(y, x, sizeof(double) * (size_t)n);
    for (int k = 0; k < stages; k++) {
        double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
        for (int i = 0; i < n; i++) {
            const double x0 = y[i];
            const double y0 = s[k].b0 * x0 + s[k].b1 * x1 + s[k].b2 * x2 - s[k].a1 * y1 - s[k].a2 * y2;
            x2 = x1; x1 = x0; y2 = y1; y1 = y0;
            y[i] = y0;
        }
    }
}

static double snr_db(const q31_t *y, const double *r, int from, int n)
{
    double ps = 0, pe = 0;
    for (int i = from; i < n; i++) {
        const double e = q31_to_double(y[i]) - r[i];
        ps += r[i] * r[i];
        pe += e * e;
    }
    return 10.0 * log10(ps / pe);
}

static void make_sine(double f, double amp)
{
    for (int i = 0; i < N; i++) {
        const double v = amp * sin(2 * M_PI * f * i / FS);
        in31[i] = q31_from_double(v);
    }
}

static void coef_quantisation(void)
{
    qfx_biquad_coef_t c;
    CHECK_EQ(qfx_biquad_coef_from_double(&c, 0.25, 0.5, 0.25, -1.9, 0.95), 0);
    CHECK_EQ(c.shift, 1);
    CHECK_EQ(c.b1, 1 << 29);
    CHECK_EQ(qfx_biquad_coef_from_double(&c, 0.1, 0.2, 0.1, -0.5, 0.2), 0);
    CHECK_EQ(c.shift, 0);
    CHECK_EQ(qfx_biquad_coef_from_double(&c, 3.0, 0, 0, 0, 0), 0);
    CHECK_EQ(c.shift, 2);
    CHECK_EQ(qfx_biquad_coef_from_double(&c, 200.0, 0, 0, 0, 0), -1);
}

static void butterworth_response(void)
{
    qfx_sos_t s[2];
    CHECK_EQ(qfx_design_butter_lowpass(FS, 1000, 4, s), 2);
    CHECK_NEAR(qfx_sos_mag_db(s, 2, FS, 0.1), 0.0, 1e-6);
    CHECK_NEAR(qfx_sos_mag_db(s, 2, FS, 1000), -3.0103, 0.01);
    /* 4th order: -24 dB/octave far above cutoff (bilinear warping adds a bit). */
    CHECK(qfx_sos_mag_db(s, 2, FS, 8000) < -70.0);
    CHECK_EQ(qfx_design_butter_lowpass(FS, 1000, 3, s), 0);
}

static int run_q31(const qfx_sos_t *s, int stages, uint8_t mode, size_t block)
{
    qfx_biquad_coef_t c[4];
    qfx_biquad_state_t st[4];
    qfx_biquad_t f;
    for (int k = 0; k < stages; k++)
        if (qfx_biquad_coef_from_double(&c[k], s[k].b0, s[k].b1, s[k].b2, s[k].a1, s[k].a2))
            return -1;
    qfx_biquad_init(&f, c, st, (uint8_t)stages, mode);
    for (size_t i = 0; i < N; i += block) {
        const size_t n = (N - i < block) ? N - i : block;
        qfx_biquad_process(&f, &in31[i], &out31[i], n);
    }
    return 0;
}

static void biquad_matches_reference(void)
{
    qfx_sos_t s[2];
    qfx_design_butter_lowpass(FS, 2000, 4, s);
    make_sine(440, 0.5);
    double x[N];
    for (int i = 0; i < N; i++) x[i] = q31_to_double(in31[i]);
    ref_cascade(s, 2, x, ref, N);

    run_q31(s, 2, QFX_BIQUAD_ROUND, N);
    CHECK(snr_db(out31, ref, 0, N) > 150.0);
    run_q31(s, 2, QFX_BIQUAD_NOISE_SHAPING, N);
    CHECK(snr_db(out31, ref, 0, N) > 150.0);
}

static void noise_shaping_helps_low_cutoff(void)
{
    /* A 20 Hz low-pass at 48 kHz has poles within 0.3 % of z = 1, where
     * quantisation noise is amplified the most. Error feedback places a
     * zero there and should win by a wide margin. */
    qfx_sos_t s[1];
    qfx_design_butter_lowpass(FS, 20, 2, s);
    make_sine(5, 0.5);
    double x[N];
    for (int i = 0; i < N; i++) x[i] = q31_to_double(in31[i]);
    ref_cascade(s, 1, x, ref, N);

    run_q31(s, 1, QFX_BIQUAD_ROUND, N);
    const double plain = snr_db(out31, ref, 0, N);
    run_q31(s, 1, QFX_BIQUAD_NOISE_SHAPING, N);
    const double shaped = snr_db(out31, ref, 0, N);
    printf("    20 Hz LP: round %.1f dB, error feedback %.1f dB\n", plain, shaped);
    CHECK(shaped > plain + 20.0);
}

static void dc_gain_is_exact_with_noise_shaping(void)
{
    qfx_sos_t s[1];
    qfx_design_butter_lowpass(FS, 50, 2, s);
    double x[N];
    for (int i = 0; i < N; i++) {
        in31[i] = q31_from_double(0.3);
        x[i] = q31_to_double(in31[i]);
    }
    ref_cascade(s, 1, x, ref, N);
    /* With error feedback the long-run average output matches the exact DC
     * gain to a fraction of an LSB; plain rounding has a fixed bias. */
    double mean_shaped = 0, mean_round = 0;
    run_q31(s, 1, QFX_BIQUAD_NOISE_SHAPING, N);
    for (int i = N / 2; i < N; i++) mean_shaped += (double)out31[i] - ref[i] * 2147483648.0;
    run_q31(s, 1, QFX_BIQUAD_ROUND, N);
    for (int i = N / 2; i < N; i++) mean_round += (double)out31[i] - ref[i] * 2147483648.0;
    mean_shaped /= N / 2;
    mean_round /= N / 2;
    printf("    DC bias: round %.2f LSB, error feedback %.3f LSB\n", mean_round, mean_shaped);
    CHECK(fabs(mean_shaped) < 0.5);
    CHECK(fabs(mean_round) > 100.0);
}

static void block_size_invariance(void)
{
    qfx_sos_t s[2];
    qfx_design_butter_lowpass(FS, 3000, 4, s);
    srand(1);
    for (int i = 0; i < N; i++)
        in31[i] = (q31_t)(((int64_t)rand() - RAND_MAX / 2) * 2);

    const size_t blocks[] = { 1, 7, 64, 1000 };
    for (int m = 0; m < 2; m++) {
        run_q31(s, 2, (uint8_t)m, N);
        memcpy(tmp31, out31, sizeof(out31));
        for (size_t b = 0; b < sizeof(blocks) / sizeof(blocks[0]); b++) {
            run_q31(s, 2, (uint8_t)m, blocks[b]);
            CHECK(memcmp(tmp31, out31, sizeof(out31)) == 0);
        }
    }
}

static void in_place_processing(void)
{
    qfx_sos_t s[1] = { qfx_design_notch(FS, 50, 10) };
    make_sine(1000, 0.4);
    run_q31(s, 1, QFX_BIQUAD_ROUND, N);
    memcpy(tmp31, out31, sizeof(out31));

    qfx_biquad_coef_t c;
    qfx_biquad_state_t st;
    qfx_biquad_t f;
    qfx_biquad_coef_from_double(&c, s[0].b0, s[0].b1, s[0].b2, s[0].a1, s[0].a2);
    qfx_biquad_init(&f, &c, &st, 1, QFX_BIQUAD_ROUND);
    qfx_biquad_process(&f, in31, in31, N);
    CHECK(memcmp(tmp31, in31, sizeof(in31)) == 0);
}

static void notch_removes_mains_hum(void)
{
    qfx_sos_t s[1] = { qfx_design_notch(FS, 50, 2) };
    make_sine(50, 0.5);
    run_q31(s, 1, QFX_BIQUAD_NOISE_SHAPING, N);
    double peak = 0;
    for (int i = N / 2; i < N; i++)
        if (fabs(q31_to_double(out31[i])) > peak) peak = fabs(q31_to_double(out31[i]));
    CHECK(20 * log10(peak / 0.5) < -40.0);
}

static void output_saturates_instead_of_wrapping(void)
{
    qfx_sos_t s[1] = { qfx_design_peaking(FS, 1000, 1, 12.0) };
    make_sine(1000, 0.9); /* +12 dB on a 0.9 sine cannot fit */
    run_q31(s, 1, QFX_BIQUAD_ROUND, N);
    int pos = 0, neg = 0;
    for (int i = N / 2; i < N; i++) {
        pos |= out31[i] == INT32_MAX;
        neg |= out31[i] == INT32_MIN;
    }
    CHECK(pos && neg);
}

static void fir_impulse_and_blocks(void)
{
    enum { T = 13 };
    q15_t h[T];
    for (int k = 0; k < T; k++) h[k] = (q15_t)(1000 * (k + 1) * (k % 3 == 0 ? -1 : 1));
    q15_t st[QFX_FIR_STATE_LEN(T)];
    qfx_fir_t f;
    qfx_fir_init(&f, h, st, T);

    q15_t x[64] = { 0 }, y[64];
    x[0] = INT16_MAX;
    qfx_fir_process(&f, x, y, 64);
    for (int k = 0; k < T; k++)
        CHECK(abs(y[k] - h[k]) <= 1);
    for (int k = T; k < 64; k++)
        CHECK_EQ(y[k], 0);

    /* Same signal, ragged blocks, in place. */
    q15_t sig[300], a[300], b[300];
    srand(7);
    for (int i = 0; i < 300; i++) sig[i] = (q15_t)(rand() % 65536 - 32768);
    qfx_fir_reset(&f);
    qfx_fir_process(&f, sig, a, 300);
    qfx_fir_reset(&f);
    memcpy(b, sig, sizeof(b));
    for (int i = 0, n = 1; i < 300; i += n, n = n % 17 + 1)
        qfx_fir_process(&f, &b[i], &b[i], (size_t)((300 - i < n) ? 300 - i : n));
    CHECK(memcmp(a, b, sizeof(a)) == 0);
}

static void cic_gain_and_wraparound(void)
{
    qfx_cic_t c;
    CHECK_EQ(qfx_cic_init(&c, 0, 8), -1);
    CHECK_EQ(qfx_cic_init(&c, 4, 16), 0);
    CHECK_EQ(qfx_cic_gain_bits(&c), 16);

    /* Full-scale DC: integrators wrap many times, output must still be exact. */
    int16_t x[4096];
    int32_t y[4096 / 16];
    for (int i = 0; i < 4096; i++) x[i] = INT16_MIN;
    CHECK_EQ(qfx_cic_process(&c, x, 4096, y), 256);
    CHECK_EQ(y[255], (int32_t)INT16_MIN * 65536);

    /* Ragged input blocks give the same stream. */
    int32_t z[256];
    size_t got = 0;
    qfx_cic_init(&c, 4, 16);
    for (int i = 0, n = 3; i < 4096; i += n, n = n * 7 % 41 + 1)
        got += qfx_cic_process(&c, &x[i], (size_t)((4096 - i < n) ? 4096 - i : n), &z[got]);
    CHECK_EQ(got, 256);
    CHECK(memcmp(y, z, sizeof(z)) == 0);
}

int main(void)
{
    RUN(coef_quantisation);
    RUN(butterworth_response);
    RUN(biquad_matches_reference);
    RUN(noise_shaping_helps_low_cutoff);
    RUN(dc_gain_is_exact_with_noise_shaping);
    RUN(block_size_invariance);
    RUN(in_place_processing);
    RUN(notch_removes_mains_hum);
    RUN(output_saturates_instead_of_wrapping);
    RUN(fir_impulse_and_blocks);
    RUN(cic_gain_and_wraparound);
    return CHECK_DONE();
}
