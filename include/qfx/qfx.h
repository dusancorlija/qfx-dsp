/*
 * qfx.h - Q15 / Q31 fixed-point primitives.
 *
 * All functions are static inline, branch-light and free of undefined
 * behaviour: intermediate results are widened before they can overflow,
 * and every narrowing step saturates.
 *
 * Right shifts of negative values are assumed to be arithmetic. This is
 * implementation-defined in C99 but true for GCC, Clang, IAR and ARMCC on
 * every target this library is meant for; tests/test_qfx.c checks it.
 */
#ifndef QFX_H
#define QFX_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int16_t q15_t; /* value = raw / 2^15, range [-1, 1) */
typedef int32_t q31_t; /* value = raw / 2^31, range [-1, 1) */

#define Q15_ONE_RAW 32768L
#define Q31_ONE_RAW 2147483648.0

/* ── saturation ─────────────────────────────────────────────────────────── */

static inline q15_t qfx_sat16(int32_t x)
{
    if (x > INT16_MAX) return INT16_MAX;
    if (x < INT16_MIN) return INT16_MIN;
    return (q15_t)x;
}

static inline q31_t qfx_sat32(int64_t x)
{
    if (x > INT32_MAX) return INT32_MAX;
    if (x < INT32_MIN) return INT32_MIN;
    return (q31_t)x;
}

/* ── conversion (host side, for tests and coefficient design) ───────────── */

static inline q15_t q15_from_double(double v)
{
    double r = v * 32768.0;
    r += (r >= 0.0) ? 0.5 : -0.5;
    if (r >= 32767.0) return INT16_MAX;
    if (r <= -32768.0) return INT16_MIN;
    return (q15_t)r;
}

static inline q31_t q31_from_double(double v)
{
    double r = v * Q31_ONE_RAW;
    r += (r >= 0.0) ? 0.5 : -0.5;
    if (r >= 2147483647.0) return INT32_MAX;
    if (r <= -2147483648.0) return INT32_MIN;
    return (q31_t)r;
}

static inline double q15_to_double(q15_t v) { return (double)v / 32768.0; }
static inline double q31_to_double(q31_t v) { return (double)v / Q31_ONE_RAW; }

/* ── arithmetic ─────────────────────────────────────────────────────────── */

static inline q15_t q15_add(q15_t a, q15_t b) { return qfx_sat16((int32_t)a + b); }
static inline q15_t q15_sub(q15_t a, q15_t b) { return qfx_sat16((int32_t)a - b); }
static inline q31_t q31_add(q31_t a, q31_t b) { return qfx_sat32((int64_t)a + b); }
static inline q31_t q31_sub(q31_t a, q31_t b) { return qfx_sat32((int64_t)a - b); }

/* Rounded multiply. (-1) * (-1) saturates to just below +1. */
static inline q15_t q15_mul(q15_t a, q15_t b)
{
    int32_t p = (int32_t)a * b;           /* Q30 */
    return qfx_sat16((p + ((int32_t)1 << 14)) >> 15);
}

static inline q31_t q31_mul(q31_t a, q31_t b)
{
    int64_t p = (int64_t)a * b;           /* Q62 */
    return qfx_sat32((p + (1LL << 30)) >> 31);
}

/* Q31 * Q15 -> Q31, the common "gain on a high-resolution signal" case. */
static inline q31_t q31_mul_q15(q31_t a, q15_t g)
{
    int64_t p = (int64_t)a * g;           /* Q46 */
    return qfx_sat32((p + (1LL << 14)) >> 15);
}

static inline q15_t q15_abs(q15_t a) { return (a == INT16_MIN) ? INT16_MAX : (q15_t)(a < 0 ? -a : a); }
static inline q31_t q31_abs(q31_t a) { return (a == INT32_MIN) ? INT32_MAX : (a < 0 ? -a : a); }

#ifdef __cplusplus
}
#endif

#endif /* QFX_H */
