/*
 * cic.h - CIC (Hogenauer) decimator, e.g. for PDM microphones and
 * sigma-delta ADC bit streams.
 *
 * Integrators run at the input rate, combs at the output rate. All state is
 * uint32_t: the integrators are allowed to wrap, and because two's-complement
 * arithmetic is modular the comb differences come out right as long as the
 * final result fits in 32 bits:
 *
 *     input_bits + order * log2(decimation) <= 32
 *
 * The DC gain is decimation^order; qfx_cic_gain_bits() returns its log2 for
 * power-of-two decimation so the output can be normalised with a shift.
 */
#ifndef QFX_CIC_H
#define QFX_CIC_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QFX_CIC_MAX_ORDER 6

typedef struct {
    uint32_t integ[QFX_CIC_MAX_ORDER];
    uint32_t comb[QFX_CIC_MAX_ORDER];
    uint16_t decimation;
    uint16_t phase;
    uint8_t order;
} qfx_cic_t;

/* Returns 0 on success, -1 on invalid parameters. */
int qfx_cic_init(qfx_cic_t *c, uint8_t order, uint16_t decimation);

/* Consumes n input samples, writes floor((phase + n) / decimation) outputs.
 * Returns the number of outputs written. */
size_t qfx_cic_process(qfx_cic_t *c, const int16_t *in, size_t n, int32_t *out);

int qfx_cic_gain_bits(const qfx_cic_t *c);

#ifdef __cplusplus
}
#endif

#endif /* QFX_CIC_H */
