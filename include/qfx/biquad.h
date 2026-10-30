/*
 * biquad.h - Q31 Direct Form I biquad cascade with optional first-order
 * error feedback (noise shaping).
 *
 *   y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2]
 *
 * Coefficients are stored as round(c * 2^(31 - shift)), so a section with
 * shift = 1 can represent coefficients in [-2, 2), which covers every stable
 * second-order denominator. Peaking / shelving sections with gain may need
 * shift = 2.
 *
 * Why DF1 and not DF2T: in fixed point, DF1 has a single quantisation point
 * per section and its states are the (bounded) input and output signals, so
 * internal overflow cannot happen. DF2T saves two state words but its states
 * can grow far beyond the signal range for high-Q, low-frequency sections.
 */
#ifndef QFX_BIQUAD_H
#define QFX_BIQUAD_H

#include <stddef.h>
#include <stdint.h>
#include "qfx.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    q31_t b0, b1, b2, a1, a2;
    uint8_t shift; /* 0..7 */
} qfx_biquad_coef_t;

typedef struct {
    q31_t x1, x2, y1, y2;
    int64_t err; /* quantisation error carried to the next sample */
} qfx_biquad_state_t;

enum {
    QFX_BIQUAD_ROUND = 0,         /* round to nearest at each section output */
    QFX_BIQUAD_NOISE_SHAPING = 1, /* first-order error feedback: zero at DC */
};

typedef struct {
    const qfx_biquad_coef_t *coef;
    qfx_biquad_state_t *state;
    uint8_t stages;
    uint8_t mode;
} qfx_biquad_t;

void qfx_biquad_init(qfx_biquad_t *f, const qfx_biquad_coef_t *coef,
                     qfx_biquad_state_t *state, uint8_t stages, uint8_t mode);

void qfx_biquad_reset(qfx_biquad_t *f);

/* Filters n samples. in and out may be the same buffer. Processing a signal
 * in blocks of any size gives bit-identical output to processing it at once. */
void qfx_biquad_process(qfx_biquad_t *f, const q31_t *in, q31_t *out, size_t n);

/* Quantise double-precision coefficients (a0 already normalised to 1).
 * Chooses the smallest shift that fits. Returns 0 on success, -1 if a
 * coefficient needs |c| >= 128. */
int qfx_biquad_coef_from_double(qfx_biquad_coef_t *c, double b0, double b1,
                                double b2, double a1, double a2);

#ifdef __cplusplus
}
#endif

#endif /* QFX_BIQUAD_H */
