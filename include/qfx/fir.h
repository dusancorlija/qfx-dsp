/*
 * fir.h - Q15 FIR filter with a mirrored delay line.
 *
 * The delay line holds every sample twice (state[i] and state[i + ntaps]),
 * so the newest ntaps samples are always contiguous in memory. The inner
 * loop is then a plain dot product with no modulo and no wrap check, which
 * compilers turn into SMLAD / SMLALD on Cortex-M4/M7.
 */
#ifndef QFX_FIR_H
#define QFX_FIR_H

#include <stddef.h>
#include <stdint.h>
#include "qfx.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const q15_t *taps; /* h[0] .. h[ntaps-1] */
    q15_t *state;      /* 2 * ntaps samples, caller-provided */
    uint16_t ntaps;
    uint16_t idx;
} qfx_fir_t;

#define QFX_FIR_STATE_LEN(ntaps) (2u * (ntaps))

void qfx_fir_init(qfx_fir_t *f, const q15_t *taps, q15_t *state, uint16_t ntaps);
void qfx_fir_reset(qfx_fir_t *f);

/* in and out may be the same buffer. Block-size invariant. */
void qfx_fir_process(qfx_fir_t *f, const q15_t *in, q15_t *out, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* QFX_FIR_H */
