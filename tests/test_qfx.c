#include "check.h"
#include "qfx/qfx.h"

static void arithmetic_shift_is_arithmetic(void)
{
    volatile int32_t a = -5;
    volatile int64_t b = -5;
    CHECK_EQ(a >> 1, -3);
    CHECK_EQ(b >> 1, -3);
}

static void saturation(void)
{
    CHECK_EQ(q15_add(INT16_MAX, 1), INT16_MAX);
    CHECK_EQ(q15_sub(INT16_MIN, 1), INT16_MIN);
    CHECK_EQ(q31_add(INT32_MAX, INT32_MAX), INT32_MAX);
    CHECK_EQ(q31_sub(INT32_MIN, 1), INT32_MIN);
    CHECK_EQ(q15_abs(INT16_MIN), INT16_MAX);
    CHECK_EQ(q31_abs(INT32_MIN), INT32_MAX);
}

static void multiply(void)
{
    /* -1 * -1 is the one product that does not fit. */
    CHECK_EQ(q15_mul(INT16_MIN, INT16_MIN), INT16_MAX);
    CHECK_EQ(q31_mul(INT32_MIN, INT32_MIN), INT32_MAX);

    CHECK_EQ(q15_mul(q15_from_double(0.5), q15_from_double(0.5)), q15_from_double(0.25));
    CHECK_EQ(q31_mul(q31_from_double(-0.5), q31_from_double(0.5)), q31_from_double(-0.25));
    CHECK_EQ(q31_mul_q15(q31_from_double(0.75), q15_from_double(0.5)), q31_from_double(0.375));

    /* Rounding is to nearest: 1 LSB * 0.5 rounds up, not down to zero. */
    CHECK_EQ(q15_mul(1, q15_from_double(0.5)), 1);
    CHECK_EQ(q15_mul(1, q15_from_double(0.25)), 0);
}

static void conversion(void)
{
    CHECK_EQ(q15_from_double(1.0), INT16_MAX);
    CHECK_EQ(q15_from_double(-1.0), INT16_MIN);
    CHECK_EQ(q31_from_double(2.0), INT32_MAX);
    CHECK_NEAR(q31_to_double(q31_from_double(0.123456789)), 0.123456789, 1e-9);
    CHECK_NEAR(q15_to_double(q15_from_double(-0.3)), -0.3, 1.0 / 32768);
}

int main(void)
{
    RUN(arithmetic_shift_is_arithmetic);
    RUN(saturation);
    RUN(multiply);
    RUN(conversion);
    return CHECK_DONE();
}
