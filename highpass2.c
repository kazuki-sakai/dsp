#include <stdio.h>
#include <errno.h>
#include "dsplib.h"

int main(int argc, char *argv[])
{
    if( argc != 4 ){
        fprintf(stderr, "Usage:\n $ %s file cutoff df\n", argv[0]);
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

    double fcut, df;
    char  *endp;
    fcut = strtod(argv[2], &endp);
    df   = strtod(argv[3], &endp);
//  fprintf(stderr, "fcut = %f\n", fcut);

    double wp[2], ws[2];
    double Gp  =  1.0;
    double As  =  6.0;
    wp[0] = (fcut + df)/nyq;
    ws[0] = (fcut - df)/nyq;
    wp[1] = -1;
    ws[1] = -1;

    dsp_sos *sos;
    sos = iir_design_sos(wp, ws, Gp, As, IIR_BUTTER);

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

    printf("# pass band edge = %f [Hz]\n", wp[0]*nyq);
    printf("# stop band edge = %f [Hz]\n", ws[0]*nyq);
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
