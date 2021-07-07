#include <stdio.h>
#include <errno.h>
#include "dsplib.h"

int main(int argc, char *argv[])
{
    /*
    dsp_zpk *zpk0;
    zpk0 = butter_zpk(5);

    freq_res_s  *fres;
    fres = freq_responce_s(zpk0, NULL, 100);

    for(uint32_t i = 0; i < fres->len; i++){
        printf("%e %e %e\n", fres->freq[i], fres->amp[i], fres->phase[i]);
    }
    return 0;
    
    dsp_butter_arg *arg;
    arg = butter_order(0.2, 0.3, 1.0, 15.0);
    printf("N = %u, Wc = %e\n", arg->order, arg->Wc);

    dsp_zpk *zpk;
    zpk = butter_zpk(5);
    //zpk = butter_zpk(arg->order);
    printf("gain = %e\n", zpk->gain);
    printf("n_zero = %u\n", zpk->n_zero);
    printf("n_pole = %u\n", zpk->n_pole);
    for(uint16_t i = 0; i < zpk->n_zero; i++){
        printf("zero[%u] = %e + I %e\n", i, creal(zpk->zero[i]), cimag(zpk->zero[i]));
    }
    for(uint16_t i = 0; i < zpk->n_pole; i++){
        printf("pole[%u] = %e + I %e\n", i, creal(zpk->pole[i]), cimag(zpk->pole[i]));
    }

    dsp_tf  *tf_a;
    tf_a = zpk2tf(zpk);

    for(uint16_t i = 0; i < tf_a->n_a; i++){
        printf("a[%u] = %e\n", i, tf_a->a[i]);
    }
    for(uint16_t i = 0; i < tf_a->n_b; i++){
        printf("b[%u] = %e\n", i, tf_a->b[i]);
    }
    return 0;

    doubleC *c;
    c = poly(zpk->pole, zpk->n_pole);
    for(uint16_t i = 0; i < zpk->n_pole + 1; i++){
        printf("c[%u] = %e + I %e\n", i, creal(c[i]), cimag(c[i]));
    }

    dsp_zpk *zpk_z;
    zpk_z = zpk_bilinear(zpk);
    printf("gain = %e\n", zpk_z->gain);
    printf("n_zero = %u\n", zpk_z->n_zero);
    printf("n_pole = %u\n", zpk_z->n_pole);
    for(uint16_t i = 0; i < zpk_z->n_zero; i++){
        printf("zero[%u] = %e + I %e\n", i, creal(zpk_z->zero[i]), cimag(zpk_z->zero[i]));
    }
    for(uint16_t i = 0; i < zpk->n_pole; i++){
        printf("pole[%u] = %e + I %e\n", i, creal(zpk_z->pole[i]), cimag(zpk_z->pole[i]));
    }


    dsp_tf *tf0;
    tf0 = zpk2tf(zpk_z);

    for(uint16_t i = 0; i < tf0->n_a; i++){
        printf("a[%u] = %e\n", i, tf0->a[i]);
    }
    for(uint16_t i = 0; i < tf0->n_b; i++){
        printf("b[%u] = %e\n", i, tf0->b[i]);
    }

    return 0;

    doubleC a[] = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    doubleC *b;

    b = poly(a, 6);
    for(uint16_t i=0; i<7; i++){
        printf("b[0]=%f + i %f\n", creal(b[i]), cimag(b[i]));
    }
    return 0;
    */

    if( argc != 3 ){
        fprintf(stderr, "Usage:\n $ %s file cutoff\n", argv[0]);
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

    double fcut;
    char  *endp;
    fcut = strtod(argv[2], &endp);
//  fprintf(stderr, "fcut = %f\n", fcut);

    double wp[2], ws[2];
    double df  = 10.0;
    double Gp  =  1.0;
    double As  = 15.0;
    wp[0] = (fcut - df)/nyq;
    ws[0] = (fcut + df)/nyq;
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
