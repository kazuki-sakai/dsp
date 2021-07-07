#pragma once

#include <errno.h>
#include <string.h>
#include "iir.h"

// #define DEBUG 1

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

#ifndef M_PIl
#define M_PIl   3.1415926535897932384626433832795029L
#endif

#define TH_COMP2REAL    1e-10
#define CPLXPAIR_EPS    1e-10

#define CALLOC_RET(x, s, t, ret)   \
                if( ((x) = (t*)calloc((s), sizeof(t))) == NULL ){  \
                       fprintf(stderr, "calloc[function=%s,line=%d]: %s\n", __func__, __LINE__, strerror(errno));  \
                       return (ret);  \
                }

#define CALLOC_GOTO(x, s, t, label)   \
                if( ((x) = (t*)calloc((s), sizeof(t))) == NULL ){  \
                       fprintf(stderr, "calloc[function=%s,line=%d]: %s\n", __func__, __LINE__, strerror(errno));  \
                       goto label;  \
                }

#define CALLOC_ABORT(x, s, t, ret)   \
                if( ((x) = (t*)calloc((s), sizeof(t))) == NULL ){  \
                       fprintf(stderr, "calloc[function=%s,line=%d]: %s\n", __func__, __LINE__, strerror(errno));  \
                       exit(EXIT_FAILURE); \
                }

#define FOPEN_RET(fp, fname, mode, ret)  \
                if( ((fp) = fopen((fname), (mode))) == NULL ){  \
                       fprintf(stderr, "fopen(%s) [function=%s,line=%d]: %s\n", (fname), __func__, __LINE__, strerror(errno));  \
                       return (ret); \
                }

#define FOPEN_GOTO(fp, fname, mode, label)  \
                if( ((fp) = fopen((fname), (mode))) == NULL ){  \
                       fprintf(stderr, "fopen(%s) [function=%s,line=%d]: %s\n", (fname), __func__, __LINE__, strerror(errno));  \
                       goto label;  \
                }

#define FOPEN_ABORT(fp, fname, mode, label)  \
                if( ((fp) = fopen((fname), (mode))) == NULL ){  \
                       fprintf(stderr, "fopen(%s) [function=%s,line=%d]: %s\n", (fname), __func__, __LINE__, strerror(errno));  \
                       exit(EXIT_FAILURE);  \
                }

#define DEBUGING  fprintf(stderr, "debug: file=%s, function=%s, line=%d\n", __FILE__, __func__, __LINE__);

/* Frequency response for analog filter */
typedef struct _freq_res_s {
    uint32_t  len;     /* length of the data */
    double   *freq;    /* frequency [Hz] */
    double   *amp;     /* amplitude */
    double   *phase;   /* phase */
} freq_res_s;

/* Frequency responce for digital filter */
typedef struct _freq_res_z {
    uint32_t  len;     /* length of the data */
    double   *w;       /* angular frequency [rad] */
    double   *amp;     /* amplitude */
    double   *phase;   /* phase */
} freq_res_z;

/** prototype **/

/**
 * Objective function of minimize filter order of butterworth band stop filter
 */
double bandstop_obj_butt (
    double  Wp,
    uint8_t index,
    double  passb[],
    double  stopb[],
    double  Gp,
    double  As
);

/**
 * To get variable to minimize the objective function ord_func
 */
double order_minimize (
    double (*ord_func)(double, uint8_t, double*, double*, double, double),
    uint8_t index,
    double  ini,
    double  fin,
    double  inc,
    double  passb[],
    double  stopb[],
    double  Gp,
    double  As
);

/**
 * Calculate nCr
 */
uint32_t nCr (
    uint32_t n,
    uint32_t r
);

/**
 * Make pairs of complex conjugates
     sorted in accordance with near unit circle
 */
bool cplxpair (
    doubleC  *arr,
    uint32_t  len
);

/**
 * Compare two complex values,
     increasing of real part
 */
int creal_incr (
    const void *a,
    const void *b
);

/**
 * Compare two complex values,
     in accordance with the distance to the unit circle
 */
int near_ucircle (
    const void *a,
    const void *b
);

/**
 * Compare two complex value,
    in accordance with the distance to the origin
 */
int near_origin (
    const void *a,
    const void *b
);

/**
 * Read file $fname for dsp
 * 1st column be stored in t
 * 2nd column be stored in x
 * Return number of data row
 */
uint32_t dsp_readfile (
    char    *fname,
    double **t,
    double **x
);

/**
 * Calculate nyquist frequency of given time series
 */
double dsp_nyqfreq (
    double   t[],
    uint32_t dlen
);

/**
 * Calculate frequency responce of analog filter
     of given [zero pole gain], or transfer function
   if both are not NULL, it uses zpk.
 */
freq_res_s *freq_responce_s (
    dsp_zpk  *zpk,    /* zeros, poles and gain */
    dsp_tf   *tf,     /* transfer function */
    double    fs      /* sampling frequency [Hz] */
);

/**
 * Calculate frequency responce of digital filter
     of given [zero pole gain], or transfer function
   if both are not NULL, it uses zpk.
 */
freq_res_z *freq_responce_z (
    dsp_zpk  *zpk,    /* zeros, poles and gain */
    dsp_tf   *tf,     /* transfer function */
    double    fs      /* sampling frequency [Hz] */
);

/**
 * Print frequency responce
 */
void print_fres_s (
    freq_res_s  *fres
);

void print_fres_z (
    freq_res_z  *fres
);

/**
 * apply linear filter to x
 */
bool linear_filter_tf (
    dsp_tf   *tf,
    uint32_t  ndata,
    double    x[],
    double    y[]
);

bool linear_filter_sos (
    dsp_sos  *sos,
    uint32_t  ndata,
    double    x[],
    double    y[]
);

/**
 * apply zero-phase filter to x
 */
bool filtfilt_tf (
    dsp_tf  *tf,
    uint32_t ndata,
    double   x[],
    double   y[]
);

bool filtfilt_sos (
    dsp_sos *sos,
    uint32_t ndata,
    double   x[],
    double   y[]
);
