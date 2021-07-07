#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include "dsplib.h"

double bandstop_obj_butt (
    double  Wp,
    uint8_t pb_index,
    double  passb[],
    double  stopb[],
    double  Gp,
    double  As
)
{
    double passband[2];

    passband[0] = ((pb_index == 0)? Wp: passb[0]);
    passband[1] = ((pb_index != 0)? Wp: passb[1]);

    double nat_t[2];
    nat_t[0] = ( stopb[0]*(passband[0] - passband[1]) )/( stopb[0]*stopb[0] - passband[0]*passband[1] );
    nat_t[1] = ( stopb[1]*(passband[0] - passband[1]) )/( stopb[1]*stopb[1] - passband[0]*passband[1] );

    double nat;
    if( fabs(nat_t[0]) <= fabs(nat_t[1]) ){
        nat = fabs(nat_t[0]);
    } else {
        nat = fabs(nat_t[1]);
    }

    double GPASS = pow(10.0, fabs(Gp)*0.1);
    double GSTOP = pow(10.0, fabs(As)*0.1);

    double n;
    n = log10( (GSTOP - 1.0)/(GPASS - 1.0) )/( 2.0*log10(nat) );

    return n;
}

double order_minimize (
    double (*ord_func)(double, uint8_t, double*, double*, double, double),
    uint8_t pb_index,
    double  ini,
    double  fin,
    double  inc,
    double  passb[],
    double  stopb[],
    double  Gp,
    double  As
)
{
    double bound;
    double min;

    min   = ord_func(ini, pb_index, passb, stopb, Gp, As);
    bound = ini;

    for(double x = ini + inc; x <= fin; x += inc){
        double y = ord_func(x, pb_index, passb, stopb, Gp, As);

        if( y < min ){
            min   = y;
            bound = x;
        }
    }

    return bound;
}

uint32_t nCr (
    uint32_t n,
    uint32_t r
)
{
    if( r == 0 ){
        return 1U;
    } else {
        uint32_t m = n * nCr(n-1, r-1);
        return m/r;
    }
}

bool cplxpair (
    doubleC  *arr,
    uint32_t  len
)
{
    doubleC  pos[len], neg[len], real[len];
    if( len%2 != 0 )    return false;

    uint32_t  i, j, k, l;
    uint32_t  p_len = 0, n_len = 0, r_len = 0;
    for(i = j = k = l = 0; i < len; i++){
        if( cimag(arr[i]) > 0 ){
            pos[p_len] = arr[i];
            p_len++;
        } else if( cimag(arr[i]) == 0 ){
            real[r_len] = arr[i];
            r_len++;
        } else {
            neg[n_len] = arr[i];
            n_len++;
        }
    }

    if( p_len != n_len )    return false;

    qsort(pos,  p_len, sizeof(doubleC), near_ucircle);
    qsort(neg,  n_len, sizeof(doubleC), near_ucircle);
    qsort(real, r_len, sizeof(doubleC), near_ucircle);

    uint8_t  *picked;
    CALLOC_RET(picked, p_len, uint8_t, false);

    doubleC  arr2[len];
    uint8_t  flg;
    for(i = 0, k = 0; i < n_len; i++){
        double r = cabs(neg[i]);
        
        for(j = 0, flg = 0; j < p_len; j++){
            if( picked[j] ) continue;
            if( ( r - cabs(pos[j]) ) <= CPLXPAIR_EPS ){
#ifdef DEBUG
                fprintf(stderr, "(%u, %u) = (%e + I*%e, %e + I*%e)\n", i,j, creal(neg[i]), cimag(neg[i]), creal(pos[j]), cimag(pos[j]));
#endif
                arr2[2*k]   = neg[i];
                arr2[2*k+1] = pos[j];
                picked[j] = 1;
                flg       = 1;
                k++;
                break;
            }
        }
        if( !flg ){
            fprintf(stderr, "not flg: %e + I*%e\n", creal(neg[i]), cimag(neg[i]));
            goto FALSE;
        }
    }

    for(i = 0; i < r_len; i++){
        arr2[2*k + i] = real[i];
    }

    for(i = 0; i < len; i++){
        arr[i] = arr2[i];
#ifdef DEBUG
        fprintf(stderr, "arr[%u] = %e + I*%e\n", i, creal(arr[i]), cimag(arr[i]));
#endif
    }

    return true;

FALSE:
    free(picked);
    return false;
}

int creal_incr (
    const void *a,
    const void *b
)
{
    double ra, rb;

    ra = creal(*((const doubleC*)a));
    rb = creal(*((const doubleC*)b));

    if( ra == rb )  return 0;
    if( ra >  rb )  return 1;
    return -1;
}

int near_origin (
    const void *a,
    const void *b
)
{
    double ra, rb;
    ra = cabs(*((const doubleC*)a));
    rb = cabs(*((const doubleC*)b));

    if( ra == rb )  return 0;
    if( ra >  rb )  return 1;
    return -1;
}

int near_ucircle (
    const void *a,
    const void *b
)
{
    double ra, rb;
    ra = fabs(1.0 - cabs(*((const doubleC*)a)));
    rb = fabs(1.0 - cabs(*((const doubleC*)b)));

    if( ra == rb ){
        double ar, br;
        ar = creal(*((const doubleC*)a));
        br = creal(*((const doubleC*)b));
        if( ar == br )  return 0;
        if( ar >  br )  return 1;
        return -1;
    }
    if( ra >  rb )  return 1;
    return -1;
}

uint32_t dsp_readfile (
    char *fname,
    double  **t,
    double  **x
)
{
    FILE *fp;
    FOPEN_RET(fp, fname, "r", 0);

    uint32_t  dlen = 0;
    char     *line = NULL;
    size_t    len  = 0;

    while( getline(&line, &len, fp) != -1 ){
        if( line[0] == '#' )    continue;
        dlen++;
    }

    CALLOC_GOTO(*t, dlen, double, FALSE1);
    CALLOC_GOTO(*x, dlen, double, FALSE2);

    fseek(fp, 0, SEEK_SET);

    for(uint32_t k = 0; k < dlen; k++){
        do{
            if( getline(&line, &len, fp) == -1 ){
                goto FALSE3;
            }
        }while( line[0] == '#' );

        double t_in, x_in;
        sscanf(line, "%lf %lf", &t_in, &x_in);
        (*t)[k] = t_in;
        (*x)[k] = x_in;
    }

    return dlen;

FALSE3:
    free(*x);
FALSE2:
    free(*t);
FALSE1:
    fclose(fp);
    return 0;
}

double dsp_nyqfreq (
    double t[],
    uint32_t dlen
)
{
    double duration, nyq;
    
    duration = t[dlen-1] - t[0];
    nyq      = (dlen - 1.0)/(2.0*duration);

    return nyq;
}

freq_res_s *freq_responce_s (
    dsp_zpk  *zpk,
    dsp_tf   *tf,
    double    fs
)
{
    if(( zpk == NULL )&&( tf == NULL )){
        fprintf(stderr, "Error: either zpk or tf must be specified.\n");
        return NULL;
    }

    freq_res_s  *fres;
    CALLOC_RET(fres, 1, freq_res_s, NULL);

    double df;
    df = 10.0/fs;
    
    fres->len = fs/2.0;
    CALLOC_GOTO(fres->freq,  fres->len,  double,  FALSE1);
    CALLOC_GOTO(fres->amp,   fres->len,  double,  FALSE2);
    CALLOC_GOTO(fres->phase, fres->len,  double,  FALSE3);

    if( zpk != NULL ){
        for(uint32_t k = 0; k < fres->len; k++){
            double f = k*df;
            double W = 2.0*M_PI*f;
            
            fres->freq[k]  = f;
            fres->amp[k]   = zpk->gain;
            fres->phase[k] = 0.0;

            for(uint16_t i = 0; i < zpk->n_zero; i++){
                fres->amp[k]   *= cabs(I*W - zpk->zero[i]);
                fres->phase[k] += carg(I*W - zpk->zero[i]);
            }
            for(uint16_t i = 0; i < zpk->n_pole; i++){
                fres->amp[k]   /= cabs(I*W - zpk->pole[i]);
                fres->phase[k] -= carg(I*W - zpk->pole[i]);
            }
        }
    }
    else {
        for(uint32_t k = 0; k < fres->len; k++){
            double f = k*df;
            double W = 2.0*M_PI*f;
            
            doubleC  numerat = 0.0;
            doubleC  denomin = 0.0;

            for(uint16_t i = 0; i < tf->n_a; i++){
                denomin += tf->a[i]*cpow(I*W, (tf->n_a - i - 1));
            }
            for(uint16_t i = 0; i < tf->n_b; i++){
                numerat += tf->b[i]*cpow(I*W, (tf->n_b - i - 1));
            }

            fres->freq[k]  = f;
            fres->amp[k]   = cabs(numerat/denomin);
            fres->phase[k] = carg(numerat/denomin);
        }
    }

    return fres;

    free(fres->phase);
FALSE3:
    free(fres->amp);
FALSE2:
    free(fres->freq);
FALSE1:
    free(fres);
    return NULL;
}

freq_res_z *freq_responce_z (
    dsp_zpk  *zpk,
    dsp_tf   *tf,
    double    fs
)
{
    if(( zpk == NULL )&&( tf == NULL )){
        fprintf(stderr, "Error: either zpk or tf must be specified.\n");
        return NULL;
    }

    freq_res_z  *fres;
    CALLOC_RET(fres, 1, freq_res_z, NULL);

    double dw;
    dw = 2.0/fs;
    
    fres->len = fs/2.0;
    CALLOC_GOTO(fres->w,     fres->len,  double,  FALSE1);
    CALLOC_GOTO(fres->amp,   fres->len,  double,  FALSE2);
    CALLOC_GOTO(fres->phase, fres->len,  double,  FALSE3);

    if( zpk != NULL ){
        for(uint32_t k = 0; k < fres->len; k++){
            double w = k*dw;
            
            fres->w[k]     = w;
            fres->amp[k]   = zpk->gain;
            fres->phase[k] = 0.0;

            for(uint16_t i = 0; i < zpk->n_zero; i++){
                fres->amp[k]   *= cabs(cexp(I*w) - zpk->zero[i]);
                fres->phase[k] += carg(cexp(I*w) - zpk->zero[i]);
            }
            for(uint16_t i = 0; i < zpk->n_pole; i++){
                fres->amp[k]   /= cabs(cexp(I*w) - zpk->pole[i]);
                fres->phase[k] -= carg(cexp(I*w) - zpk->pole[i]);
            }
        }
    }
    else {
        for(uint32_t k = 0; k < fres->len; k++){
            double w = k*dw;
            
            doubleC  numerat = 0.0;
            doubleC  denomin = 0.0;

            for(uint16_t i = 0; i < tf->n_a; i++){
                denomin += tf->a[i]*cexp(-I*w*i);
            }
            for(uint16_t i = 0; i < tf->n_b; i++){
                numerat += tf->b[i]*cexp(-I*w*i);
            }

            fres->w[k]     = w;
            fres->amp[k]   = cabs(numerat/denomin);
            fres->phase[k] = carg(numerat/denomin);
        }
    }

    return fres;

    free(fres->phase);
FALSE3:
    free(fres->amp);
FALSE2:
    free(fres->w);
FALSE1:
    free(fres);
    return NULL;
}

void print_fres_s (
    freq_res_s  *fres
)
{
    for(uint32_t k = 0; k < fres->len; k++){
        printf("%e %e %e\n", fres->freq[k], fres->amp[k], fres->phase[k]);
    }
}

void print_fres_z (
    freq_res_z  *fres
)
{
    for(uint32_t k = 0; k < fres->len; k++){
        printf("%e %e %e\n", fres->w[k], fres->amp[k], fres->phase[k]);
    }
}


bool linear_filter_tf (
    dsp_tf  *tf,
    uint32_t ndata,
    double   x[],
    double   y[]
)
{
    uint16_t n_f = ((tf->n_a >= tf->n_b)? tf->n_a: tf->n_b);
    double   z[n_f], z0[n_f];

    //fprintf(stderr, "a[0] = %e\n", tf->a[0]);
    double a[n_f], b[n_f];
    for(uint16_t m = 0; m < n_f; m++){
        a[m] = (( m < tf->n_a)? (tf->a[m])/(tf->a[0]): 0);
        b[m] = (( m < tf->n_b)? (tf->b[m])/(tf->a[0]): 0);
    }

    for(uint32_t m = 0; m < n_f; m++){
        z0[m] = 0.0;
    }
    z[n_f - 1] = 0.0;

    for(uint32_t m = 0; m < ndata; m++){
        y[m]  = b[0]*x[m] + z0[0];
        for(uint32_t n = 0; n < n_f - 1; n++){
            z[n]  = b[n+1]*x[m] + z0[n+1] - a[n+1]*y[m];
            z0[n] = z[n];
        }
    }

    return true;
}

bool linear_filter_sos (
    dsp_sos *sos,
    uint32_t ndata,
    double   x[],
    double   y[]
)
{
    double x1, x2, y1, y2, *w;

    double **a;
    a = sos->a;

    CALLOC_RET(w, ndata, double, false);

    for(uint32_t n = 0; n < ndata; n++){
        w[n] = x[n];
    }

    for(uint16_t i = 0; i < sos->L; i++){
        x1 = x2 = y1 = y2 = 0.0;
        for(uint32_t n = 0; n < ndata; n++){
            y[n] = (a[i][0]*w[n] + a[i][1]*x1 + a[i][2]*x2 - a[i][4]*y1 - a[i][5]*y2)/a[i][3];
            x2 = x1;
            x1 = w[n];
            y2 = y1;
            y1 = y[n];
            w[n] = y[n];
        }
    }

    for(uint32_t n = 0; n < ndata; n++){
        y[n] *= sos->G;
    }

    free(w);
    return true;
}

bool filtfilt_tf (
    dsp_tf  *tf,
    uint32_t ndata,
    double   x[],
    double   y[]
)
{
    double  y_t[ndata];
    double  y0[ndata], y1[ndata];

    linear_filter_tf(tf, ndata, x, y_t);

    for(uint32_t n = 0; n < ndata; n++){
        y0[n] = y_t[ndata - n - 1];
    }

    linear_filter_tf(tf, ndata, y0, y1);

    for(uint32_t n = 0; n < ndata; n++){
        y[n]  = y1[ndata - n - 1];
    }

    return true;
}

bool filtfilt_sos (
    dsp_sos  *sos,
    uint32_t  ndata,
    double    x[],
    double    y[]
)
{
    double *y_t = NULL;
    double *y0  = NULL;
    double *y1  = NULL;
    CALLOC_GOTO(y_t, ndata, double, FALSE);
    CALLOC_GOTO(y0, ndata, double, FALSE);
    CALLOC_GOTO(y1, ndata, double, FALSE);

    linear_filter_sos(sos, ndata, x, y_t);

    for(uint32_t n = 0; n < ndata; n++){
        y0[n] = y_t[ndata - n - 1];
    }

    linear_filter_sos(sos, ndata, y0, y1);

    for(uint32_t n = 0; n < ndata; n++){
        y[n] = y1[ndata - n - 1];
    }

    return true;

FALSE:
    if( y_t )   free(y_t);
    if( y0  )   free(y0);
    if( y1  )   free(y1);

    return false;
}
