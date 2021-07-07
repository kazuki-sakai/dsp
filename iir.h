/********************************/
/*            dsp.h             */
/*                   ver. 1.0   */
/*           (c) Kazuki SAKAI   */
/********************************/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

typedef double _Complex doubleC;

typedef enum {
    LOWPASS,
    HIGHPASS,
    BANDPASS,
    BANDSTOP,
} filter_type;

typedef enum {
    IIR_BUTTER,
    IIR_ELLIP,
} iir_type;

/**
 * zeloes, poles and gain of transfer function
 */
typedef struct _dsp_zpk {
    uint16_t  n_zero;
    uint16_t  n_pole;
    double    gain;
    doubleC  *zero;
    doubleC  *pole;
} dsp_zpk;

/**
 * Coefficients of numerator and denominator of transfer function
 */
typedef struct _dsp_tf {
    uint16_t  n_a;
    uint16_t  n_b;
    double   *a;
    double   *b;
} dsp_tf;

/**
 * Coefficients matrix of transfer function in the form of second order sections (Biquad)
 */
typedef struct _dsp_sos {
    uint16_t  L;
    double  **a;
    double    G;
} dsp_sos;

/**
 * Arguments of butterworth lowpass filter
 */
typedef struct _dsp_butter_arg {
    uint16_t  order;  /* filter order */
    double    wc[2];  /* cutoff frequencies [rad] */
} dsp_butter_arg;

/*******************************************/
/*   prototype declarations of functions   */
/*******************************************/

/**
 * Calculating convolution
 */
doubleC *convolve (
    doubleC  a[],
    doubleC  v[],
    uint16_t n_a,
    uint16_t n_v
);

/**
 * Calculating coefficients of the polynomial with given zeros
 */
doubleC *poly (
    doubleC  zero[],
    uint16_t n_zero
);

/**
 * Calcurating filter cutoff frequency of butterwise filter
 */
dsp_butter_arg *butter_cutoff (
    double wp[], /* passband edges */
    double ws[], /* stopband edges */
    double Gp,   /* gain at pass band edge */
    double As,   /* attenuation at stop band edge */
    uint16_t N,  /* filter order */
    filter_type ftype
);

/**
 * Calculating filter order and cutoff frequency of butterwise filter
 */
dsp_butter_arg *butter_order (
    double wp[], /* passband edges */
    double ws[], /* stopband edges */
    double Gp,   /* gain at pass band edge */
    double As,   /* attenuation at stop band edge */
    filter_type ftype
);

/**
 * Calculating zeros, poles and gain of butterwise filter with given order
 */
dsp_zpk *butter_zpk (
    uint16_t N
);

/**
 * Calculating zeros, poles and gain of butterwise with given cutoff frequency from given zpk of base filter
 */
dsp_zpk *zpklp2lp (
    dsp_zpk *zpk,  /* zpk of base filter */
    double   Wc    /* Prewarped cutoff frequency */
);

/**
 * Calculating zeros, poles and gain of butterwise with given cutoff frequency from given zpk of base filter
 */
dsp_zpk *zpklp2hp (
    dsp_zpk *zpk,  /* zpk of base filter */
    double   Wc    /* Prewarped cutoff frequency */
);

/**
 * Calculating zeros, poles and gain of butterwise with given cutoff frequency from given zpk of base filter
 */
dsp_zpk *zpklp2bp (
    dsp_zpk *zpk,  /* zpk of base filter */
    double   Wc[]  /* Prewarped cutoff frequencies */
);

/**
 * Calculating zeros, poles and gain of butterwise with given cutoff frequency from given zpk of base filter
 */
dsp_zpk *zpklp2bs (
    dsp_zpk *zpk,  /* zpk of base filter */
    double   Wc[]  /* Prewarped cutoff frequencies */
);

/**
 * Warping zpk of analog filter to digital filter by bilinear warping
 */
dsp_zpk *zpk_bilinear (
    dsp_zpk *zpk  /* zpk of analog filter */
);

/**
 * Get coefficients of transfer function with given zeors and poles
 */
dsp_tf *zpk2tf (
    dsp_zpk *zkp   /* zpk of transfer function */
);

/**
 * Convert zpk representation of transfer function into sos representation
 */
dsp_sos *zpk2sos (
    dsp_zpk  *zpk
);

/**
 * Warping tf of analog filter to digital filter by bilinear warping
 */
dsp_tf *bilinear (
    dsp_tf  *tf_s
);

/**
 * IIR Filter design
 *   Returning zpk of the transfer function
 */
dsp_zpk *iir_design_zpk (
    double   wp[],    /* pass band edges. wp[0] is lower edge, wp[1] is higher edge. */
    double   ws[],    /* stop band edges. ws[0] is lower edge, ws[1] is higher edge */
    double   gpass,   /* maximum loss in the passband */
    double   gstop,   /* minimum attenuation in the stopband */
    iir_type type     /* type of filter design */
);

dsp_sos *iir_design_sos (
    double   wp[],
    double   ws[],
    double   Gp,
    double   As,
    iir_type type
);

void free_dsp_tf (
    dsp_tf *tf
);

void free_dsp_sos (
    dsp_sos *sos
);
