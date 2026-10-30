#include "qfx/biquad.h"

#include <math.h>
#include <string.h>

void qfx_biquad_init(qfx_biquad_t *f, const qfx_biquad_coef_t *coef,
                     qfx_biquad_state_t *state, uint8_t stages, uint8_t mode)
{
    f->coef = coef;
    f->state = state;
    f->stages = stages;
    f->mode = mode;
    qfx_biquad_reset(f);
}

void qfx_biquad_reset(qfx_biquad_t *f)
{
    memset(f->state, 0, sizeof(*f->state) * f->stages);
}

/* Multiply-accumulate in uint64_t: signed overflow is undefined in C, but
 * unsigned wrap-around is not. The partial sums of five Q62 products can
 * exceed int64_t even when the final sum fits; modular arithmetic makes
 * those transient overflows harmless. */
#define MAC(acc, c, s) ((acc) += (uint64_t)((int64_t)(c) * (int64_t)(s)))
#define MSC(acc, c, s) ((acc) -= (uint64_t)((int64_t)(c) * (int64_t)(s)))

static void section_round(const qfx_biquad_coef_t *c, qfx_biquad_state_t *s,
                          const q31_t *in, q31_t *out, size_t n)
{
    const q31_t b0 = c->b0, b1 = c->b1, b2 = c->b2, a1 = c->a1, a2 = c->a2;
    const int sh = 31 - c->shift;
    const uint64_t half = (uint64_t)1 << (sh - 1);
    q31_t x1 = s->x1, x2 = s->x2, y1 = s->y1, y2 = s->y2;

    for (size_t i = 0; i < n; i++) {
        const q31_t x0 = in[i];
        uint64_t acc = half;
        MAC(acc, b0, x0);
        MAC(acc, b1, x1);
        MAC(acc, b2, x2);
        MSC(acc, a1, y1);
        MSC(acc, a2, y2);
        const q31_t y0 = qfx_sat32((int64_t)acc >> sh);
        x2 = x1; x1 = x0;
        y2 = y1; y1 = y0;
        out[i] = y0;
    }
    s->x1 = x1; s->x2 = x2; s->y1 = y1; s->y2 = y2;
}

static void section_shaped(const qfx_biquad_coef_t *c, qfx_biquad_state_t *s,
                           const q31_t *in, q31_t *out, size_t n)
{
    const q31_t b0 = c->b0, b1 = c->b1, b2 = c->b2, a1 = c->a1, a2 = c->a2;
    const int sh = 31 - c->shift;
    const int64_t unit = (int64_t)1 << sh;
    q31_t x1 = s->x1, x2 = s->x2, y1 = s->y1, y2 = s->y2;
    int64_t err = s->err;

    for (size_t i = 0; i < n; i++) {
        const q31_t x0 = in[i];
        uint64_t acc = (uint64_t)err;
        MAC(acc, b0, x0);
        MAC(acc, b1, x1);
        MAC(acc, b2, x2);
        MSC(acc, a1, y1);
        MSC(acc, a2, y2);
        const int64_t sacc = (int64_t)acc;
        const int64_t yq = sacc >> sh;       /* floor */
        err = sacc - yq * unit;              /* in [0, 2^sh): fed back next sample */
        const q31_t y0 = qfx_sat32(yq);
        x2 = x1; x1 = x0;
        y2 = y1; y1 = y0;
        out[i] = y0;
    }
    s->x1 = x1; s->x2 = x2; s->y1 = y1; s->y2 = y2;
    s->err = err;
}

void qfx_biquad_process(qfx_biquad_t *f, const q31_t *in, q31_t *out, size_t n)
{
    /* Stage-outer loop: each section's coefficients and state stay in
     * registers for the whole block instead of being reloaded per sample. */
    const q31_t *src = in;
    for (uint8_t k = 0; k < f->stages; k++) {
        if (f->mode == QFX_BIQUAD_NOISE_SHAPING)
            section_shaped(&f->coef[k], &f->state[k], src, out, n);
        else
            section_round(&f->coef[k], &f->state[k], src, out, n);
        src = out;
    }
    if (f->stages == 0 && in != out)
        memmove(out, in, n * sizeof(*out));
}

int qfx_biquad_coef_from_double(qfx_biquad_coef_t *c, double b0, double b1,
                                double b2, double a1, double a2)
{
    const double v[5] = { b0, b1, b2, a1, a2 };
    double m = 0.0;
    for (int i = 0; i < 5; i++)
        if (fabs(v[i]) > m) m = fabs(v[i]);

    uint8_t shift = 0;
    while (shift < 8 && m >= ldexp(1.0, shift) * (1.0 - 1e-12))
        shift++;
    if (shift >= 8)
        return -1;

    const double scale = ldexp(1.0, 31 - shift);
    q31_t *dst[5] = { &c->b0, &c->b1, &c->b2, &c->a1, &c->a2 };
    for (int i = 0; i < 5; i++) {
        double r = v[i] * scale;
        r += (r >= 0.0) ? 0.5 : -0.5;
        if (r > 2147483647.0) r = 2147483647.0;
        if (r < -2147483648.0) r = -2147483648.0;
        *dst[i] = (q31_t)r;
    }
    c->shift = shift;
    return 0;
}
