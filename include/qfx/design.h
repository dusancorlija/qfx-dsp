/*
 * design.h - double-precision biquad design (host side or FPU targets).
 *
 * Formulas follow the well-known RBJ "Audio EQ Cookbook". Results are
 * normalised so that a0 = 1 and are returned as { b0, b1, b2, a1, a2 }.
 */
#ifndef QFX_DESIGN_H
#define QFX_DESIGN_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double b0, b1, b2, a1, a2; } qfx_sos_t;

qfx_sos_t qfx_design_lowpass(double fs, double f0, double q);
qfx_sos_t qfx_design_highpass(double fs, double f0, double q);
qfx_sos_t qfx_design_bandpass(double fs, double f0, double q);
qfx_sos_t qfx_design_notch(double fs, double f0, double q);
qfx_sos_t qfx_design_peaking(double fs, double f0, double q, double gain_db);

/* Butterworth low-pass of even order as order/2 cascaded sections.
 * Writes order/2 sections to sos. Returns the number written, or 0 if the
 * order is odd or zero. */
int qfx_design_butter_lowpass(double fs, double f0, int order, qfx_sos_t *sos);

/* Magnitude response of a cascade at frequency f, in dB. */
double qfx_sos_mag_db(const qfx_sos_t *sos, int n, double fs, double f);

#ifdef __cplusplus
}
#endif

#endif /* QFX_DESIGN_H */
