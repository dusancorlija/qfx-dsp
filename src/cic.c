#include "qfx/cic.h"

#include <string.h>

int qfx_cic_init(qfx_cic_t *c, uint8_t order, uint16_t decimation)
{
    if (order == 0 || order > QFX_CIC_MAX_ORDER || decimation == 0)
        return -1;
    memset(c, 0, sizeof(*c));
    c->order = order;
    c->decimation = decimation;
    return 0;
}

size_t qfx_cic_process(qfx_cic_t *c, const int16_t *in, size_t n, int32_t *out)
{
    const uint8_t N = c->order;
    size_t produced = 0;

    for (size_t i = 0; i < n; i++) {
        uint32_t v = (uint32_t)(int32_t)in[i];
        for (uint8_t k = 0; k < N; k++) {
            c->integ[k] += v;
            v = c->integ[k];
        }
        if (++c->phase < c->decimation)
            continue;
        c->phase = 0;

        for (uint8_t k = 0; k < N; k++) {
            const uint32_t prev = c->comb[k];
            c->comb[k] = v;
            v -= prev;
        }
        out[produced++] = (int32_t)v;
    }
    return produced;
}

int qfx_cic_gain_bits(const qfx_cic_t *c)
{
    int bits = 0;
    uint16_t r = c->decimation;
    while (r > 1) {
        if (r & 1)
            return -1; /* not a power of two */
        r >>= 1;
        bits++;
    }
    return bits * c->order;
}
