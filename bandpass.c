#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "dsplib.h"

int main(int argc, char *argv[])
{
    if( argc != 4 ){
        fprintf(stderr, "Usage:\n $ %s file low_cutoff high_cutoff\n", argv[0]);
        return 1;
    }

    double *t, *x;
    uint32_t dlen;
    dlen = dsp_readfile(argv[1], &t, &x);
    if( dlen == 0 ){
        fprintf(stderr, "Error: irregal format data file '%s'\n", argv[1]);
        return 1;
    }

    double *y;
    CALLOC_RET(y, dlen, double, 1);

    double nyq;
    nyq = dsp_nyqfreq(t, dlen);
//  fprintf(stderr, "nyq = %f\n", nyq);

    double fcut[2];
    char  *endp;
    fcut[0] = strtod(argv[2], &endp);
    fcut[1] = strtod(argv[3], &endp);
//  fprintf(stderr, "fcut = %f\n", fcut);

    double wp[2], ws[2];
    double df  =  10.0;
    double Gp  =  1.0;
    double As  = 15.0;
    wp[0] = (fcut[0] + df)/nyq;
    ws[0] = (fcut[0] - df)/nyq;
    wp[1] = (fcut[1] - df)/nyq;
    ws[1] = (fcut[1] + df)/nyq;

    dsp_sos *sos;
//  sos = iir_design_sos(wp, ws, Gp, As, IIR_BUTTER);

    double fs = 2.0;
    double Wc[2];
    Wc[0] = 2.0*fs*tan(M_PI*fcut[0]/nyq/fs);
    Wc[1] = 2.0*fs*tan(M_PI*fcut[1]/nyq/fs);

    dsp_zpk *zpklp;
    zpklp = butter_zpk(10);

    dsp_zpk *zpkbp;
    zpkbp = zpklp2bp(zpklp, Wc);

    dsp_zpk *zpk_z;
    zpk_z = zpk_bilinear(zpkbp);

    sos = zpk2sos(zpk_z);

/*
    fprintf(stderr, "n_a:%u\nn_b:%u\n", tf->n_a, tf->n_b);
    for(uint16_t i = 0; i < tf->n_a; i++){
        fprintf(stderr, "a[%u] = %e\n", i, tf->a[i]);
    }
    for(uint16_t i = 0; i < tf->n_b; i++){
        fprintf(stderr, "b[%u] = %e\n", i, tf->b[i]);
    }
*/

    filtfilt_sos(sos, dlen, x, y);

    fprintf(stderr, "debug\n");

    printf("# pass band edge = %f : %f [Hz]\n", wp[0]*nyq, wp[1]*nyq);
    printf("# stop band edge = %f , %f [Hz]\n", ws[0]*nyq, ws[1]*nyq);
    printf("# pass band limit attenuation   = %f [dB]\n", Gp);
    printf("# stop band maximum attenuation = %f [dB]\n", As);
    for(uint32_t k = 0; k < dlen; k++){
        printf("%.10e %.10e\n", t[k], y[k]);
    }

/*
    freq_res_z  *fres;
    fres = freq_responce_z(NULL, tf, nyq*2);
    for(uint32_t i = 0; i < fres->len; i++){
        fprintf(stderr, "%e %e %e\n", fres->w[i]*nyq, fres->amp[i], fres->phase[i]);
    }
*/

    free(y);
    free_dsp_sos(sos);

    return 0;
}
