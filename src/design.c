#include "qfx/design.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static qfx_sos_t normalise(double b0, double b1, double b2,
                           double a0, double a1, double a2)
{
    qfx_sos_t s = { b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0 };
    return s;
}

qfx_sos_t qfx_design_lowpass(double fs, double f0, double q)
{
    const double w = 2.0 * M_PI * f0 / fs, cw = cos(w), al = sin(w) / (2.0 * q);
    return normalise((1 - cw) / 2, 1 - cw, (1 - cw) / 2, 1 + al, -2 * cw, 1 - al);
}

qfx_sos_t qfx_design_highpass(double fs, double f0, double q)
{
    const double w = 2.0 * M_PI * f0 / fs, cw = cos(w), al = sin(w) / (2.0 * q);
    return normalise((1 + cw) / 2, -(1 + cw), (1 + cw) / 2, 1 + al, -2 * cw, 1 - al);
}

qfx_sos_t qfx_design_bandpass(double fs, double f0, double q)
{
    const double w = 2.0 * M_PI * f0 / fs, cw = cos(w), al = sin(w) / (2.0 * q);
    return normalise(al, 0, -al, 1 + al, -2 * cw, 1 - al);
}

qfx_sos_t qfx_design_notch(double fs, double f0, double q)
{
    const double w = 2.0 * M_PI * f0 / fs, cw = cos(w), al = sin(w) / (2.0 * q);
    return normalise(1, -2 * cw, 1, 1 + al, -2 * cw, 1 - al);
}

qfx_sos_t qfx_design_peaking(double fs, double f0, double q, double gain_db)
{
    const double A = pow(10.0, gain_db / 40.0);
    const double w = 2.0 * M_PI * f0 / fs, cw = cos(w), al = sin(w) / (2.0 * q);
    return normalise(1 + al * A, -2 * cw, 1 - al * A, 1 + al / A, -2 * cw, 1 - al / A);
}

int qfx_design_butter_lowpass(double fs, double f0, int order, qfx_sos_t *sos)
{
    if (order <= 0 || (order & 1))
        return 0;
    const int n = order / 2;
    for (int k = 0; k < n; k++) {
        /* Pole pair angle of an order-N Butterworth prototype. */
        const double theta = M_PI * (2.0 * k + 1.0) / (2.0 * order);
        const double q = 1.0 / (2.0 * sin(theta));
        sos[k] = qfx_design_lowpass(fs, f0, q);
    }
    return n;
}

double qfx_sos_mag_db(const qfx_sos_t *sos, int n, double fs, double f)
{
    const double w = 2.0 * M_PI * f / fs;
    const double c1 = cos(w), s1 = sin(w), c2 = cos(2 * w), s2 = sin(2 * w);
    double mag = 1.0;
    for (int k = 0; k < n; k++) {
        const qfx_sos_t *s = &sos[k];
        const double nr = s->b0 + s->b1 * c1 + s->b2 * c2;
        const double ni = -(s->b1 * s1 + s->b2 * s2);
        const double dr = 1.0 + s->a1 * c1 + s->a2 * c2;
        const double di = -(s->a1 * s1 + s->a2 * s2);
        mag *= sqrt((nr * nr + ni * ni) / (dr * dr + di * di));
    }
    return 20.0 * log10(mag);
}
