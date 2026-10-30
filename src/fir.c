#include "qfx/fir.h"

#include <string.h>

void qfx_fir_init(qfx_fir_t *f, const q15_t *taps, q15_t *state, uint16_t ntaps)
{
    f->taps = taps;
    f->state = state;
    f->ntaps = ntaps;
    qfx_fir_reset(f);
}

void qfx_fir_reset(qfx_fir_t *f)
{
    memset(f->state, 0, QFX_FIR_STATE_LEN(f->ntaps) * sizeof(q15_t));
    f->idx = 0;
}

void qfx_fir_process(qfx_fir_t *f, const q15_t *in, q15_t *out, size_t n)
{
    const uint16_t N = f->ntaps;
    const q15_t *h = f->taps;
    q15_t *st = f->state;
    uint16_t idx = f->idx;

    for (size_t i = 0; i < n; i++) {
        idx = (idx == 0) ? (uint16_t)(N - 1) : (uint16_t)(idx - 1);
        st[idx] = st[idx + N] = in[i];

        /* st[idx + k] == x[n - k]. Four independent accumulators break the
         * dependency chain so the multiplier pipeline stays full. */
        const q15_t *x = &st[idx];
        int64_t a0 = 0, a1 = 0, a2 = 0, a3 = 0;
        uint16_t k = 0;
        for (; k + 4 <= N; k += 4) {
            a0 += (int32_t)h[k + 0] * x[k + 0];
            a1 += (int32_t)h[k + 1] * x[k + 1];
            a2 += (int32_t)h[k + 2] * x[k + 2];
            a3 += (int32_t)h[k + 3] * x[k + 3];
        }
        for (; k < N; k++)
            a0 += (int32_t)h[k] * x[k];

        const int64_t acc = a0 + a1 + a2 + a3;          /* Q30 */
        out[i] = qfx_sat16((int32_t)qfx_sat32((acc + (1 << 14)) >> 15));
    }
    f->idx = idx;
}
