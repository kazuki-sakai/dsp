#include <math.h>
#include <errno.h>
#include "dsplib.h"
#include "iir.h"

doubleC *convolve (
    doubleC a[],
    doubleC v[],
    uint16_t  n_a,
    uint16_t  n_v
)
{
    uint16_t n_b = n_a + n_v - 1;
    doubleC *b;

    CALLOC_RET(b, n_b, doubleC, NULL);

    for(uint16_t n = 0; n < n_b; n++){
        b[n] = 0;

        uint16_t m_fin = (((n + 1) <= n_a)? (n + 1): n_a);
        for(uint16_t m = 0; m < m_fin; m++){
            if( (n - m) >= n_v )    continue;
            b[n] += a[m]*v[n - m];
        }
    }

    return b;
}

doubleC *poly (
    doubleC zero[],
    uint16_t n_zero
)
{
    doubleC *a, v[2];

    CALLOC_GOTO(a, (n_zero + 1), doubleC, FALSE1);

    a[0] = 1.0;
    v[0] = 1.0;
    for(uint16_t k = 0; k < n_zero; k++){
        v[1] = -zero[k];

        doubleC *b;
        b = convolve(a, v, (k + 1), 2);
        if( b == NULL ) goto FALSE2;

        for(uint16_t i = 0; i < (2 + k); i++){
            a[i] = b[i];
        }
        free(b);
    }

    return a;

FALSE2:
    if( a ) free(a);
FALSE1:
    return NULL;
}

dsp_zpk *iir_design_zpk (
    double   wp[],
    double   ws[],
    double   gpass,
    double   gstop,
    iir_type type
)
{
    return NULL;
}


dsp_butter_arg *butter_cutoff (
    double  wp[],
    double  ws[],
    double  Gp,
    double  As,
    uint16_t N,
    filter_type  ftype
)
{
    dsp_butter_arg *butarg;
    CALLOC_RET(butarg, 1, dsp_butter_arg, NULL);

    /* Pre-warping */
    double Wp[2] = {0};
    // double Ws[2] = {0};
    Wp[0] = 2.0*tan(M_PI*wp[0]/2.0);
    // Ws[0] = 2.0*tan(M_PI*ws[0]/2.0);
    if(( ftype == BANDPASS )||( ftype == BANDSTOP )){
        Wp[1] = 2.0*tan(M_PI*wp[1]/2.0);
        // Ws[1] = 2.0*tan(M_PI*ws[1]/2.0);
    }

    // double G, A;
    double G;
    G = pow(10.0, Gp*0.1);
    // A = pow(10.0, As*0.1);

    double W0, Wc[2] = {0};
    W0 = pow((G - 1.0), -1.0/(2.0*N));

    double discr;
    switch( ftype ){
        case LOWPASS:
            Wc[0] = Wp[0]*W0;
            Wc[1] = -1;
            break;

        case HIGHPASS:
            Wc[0] = Wp[0]/W0;
            Wc[1] = -1;
            break;

        case BANDPASS:
            Wc[0] = ( W0*(Wp[1] - Wp[0])/2.0 + sqrt((W0*W0/4.0)*(Wp[1] - Wp[0])*(Wp[1] - Wp[0]) + Wp[0]*Wp[1]));
            Wc[1] = (-W0*(Wp[1] - Wp[0])/2.0 + sqrt((W0*W0/4.0)*(Wp[1] - Wp[0])*(Wp[1] - Wp[0]) + Wp[0]*Wp[1]));

            if( fabs(Wc[1]) < fabs(Wc[0]) ){
                double tmp;
                tmp   = fabs(Wc[0]);
                Wc[0] = fabs(Wc[1]);
                Wc[1] = tmp;
            }
            break;

        case BANDSTOP:
            discr = sqrt( (Wp[0] - Wp[1])*(Wp[0] - Wp[1]) + 4.0*W0*W0*Wp[0]*Wp[1] );
            Wc[0] = ( Wp[1] - Wp[0] + discr )/(2.0*W0);
            Wc[1] = ( Wp[1] - Wp[0] - discr )/(2.0*W0);

            if( fabs(Wc[1]) < fabs(Wc[0]) ){
                double tmp;
                tmp   = fabs(Wc[0]);
                Wc[0] = fabs(Wc[1]);
                Wc[1] = tmp;
            }
            break;
    }

    /* Warping */
    double wc[2];
    wc[0] = (2.0/M_PI)*atan(Wc[0]/2.0);
    wc[1] = (2.0/M_PI)*atan(Wc[1]/2.0);

    butarg->order = N;
    butarg->wc[0] = wc[0];
    butarg->wc[1] = wc[1];

    return butarg;
}


dsp_butter_arg *butter_order (
    double  wp[],
    double  ws[],
    double  Gp,
    double  As,
    filter_type  ftype
)
{

    /* Pre-warping */
    double Wp[2], Ws[2];
    Wp[0] = 2.0*tan(M_PI*wp[0]/2.0);
    Ws[0] = 2.0*tan(M_PI*ws[0]/2.0);
    if(( ftype == BANDPASS )||( ftype == BANDSTOP )){
        Wp[1] = 2.0*tan(M_PI*wp[1]/2.0);
        Ws[1] = 2.0*tan(M_PI*ws[1]/2.0);
    }

    double G, A;
    G = pow(10.0, Gp*0.1);
    A = pow(10.0, As*0.1);

    /* W at calculating N */
    double W_N;
    double W_N1, W_N2;
    double Wp0, Wp1;
    switch( ftype ){
        case LOWPASS:
            W_N = fabs(Ws[0]/Wp[0]);
            break;

        case HIGHPASS:
            W_N = fabs(Wp[0]/Ws[0]);
            break;

        case BANDPASS:
            W_N1 = ((Ws[0]*Ws[0] - Wp[0]*Wp[1]) / (Ws[0]*(Wp[0] - Wp[1])));
            W_N2 = ((Ws[1]*Ws[1] - Wp[0]*Wp[1]) / (Ws[1]*(Wp[0] - Wp[1])));
            
            if( fabs(W_N1) <= fabs(W_N2) ){
                W_N = fabs(W_N1);
            } else {
                W_N = fabs(W_N2);
            }
            break;

        case BANDSTOP:
            Wp0 = order_minimize(&bandstop_obj_butt, 0, Wp[0], Ws[0] - 1e-21, 1e-5, Wp, Ws, Gp, As);
            Wp[0] = Wp0;

            Wp1 = order_minimize(&bandstop_obj_butt, 1, Ws[1] + 1e-21, Wp[1], 1e-5, Wp, Ws, Gp, As);
            Wp[1] = Wp1;

            W_N1 = ( Ws[0]*(Wp[0] - Wp[1]) )/( Ws[0]*Ws[0] - Wp[0]*Wp[1] );
            W_N2 = ( Ws[1]*(Wp[0] - Wp[1]) )/( Ws[1]*Ws[1] - Wp[0]*Wp[1] );

            if( fabs(W_N1) <= fabs(W_N2) ){
                W_N = fabs(W_N1);
            } else {
                W_N = fabs(W_N2);
            }
            break;

        default:
            W_N = 0.0;
    }

    /* filter order */
    uint16_t N;
    N = (uint16_t)ceil( log10( (A - 1.0)/(G - 1.0) ) / (2.0*log10(W_N)) );
    N += N%2;

    dsp_butter_arg *butarg;
    butarg = butter_cutoff(wp, ws, Gp, As, N, ftype);
    if( butarg == NULL )    return NULL;

    return butarg;
}

dsp_zpk *butter_zpk (
    uint16_t N
)
{
    dsp_zpk *zpk;
    CALLOC_RET(zpk, 1, dsp_zpk, NULL);

    zpk->n_zero = 0;
    zpk->n_pole = N;
    zpk->gain   = 1.0;
    zpk->zero   = NULL;

    CALLOC_GOTO(zpk->pole, N, doubleC, FALSE1);

    for(uint16_t i = 0; i < N; i++){
        double theta = (double)(2*i + N + 1)/(2*N)*M_PI;

        zpk->pole[i] = cos(theta);
        if( fabs(sin(theta)) >= TH_COMP2REAL ){
            zpk->pole[i] += I*sin(theta);
        }
#ifdef DEBUG
        fprintf(stderr, "pole[%u] = %e + I*%e\n", i, creal(zpk->pole[i]), cimag(zpk->pole[i]));
#endif
    }

    return zpk;

FALSE1:
    free(zpk);
    return NULL;
}

dsp_zpk *zpklp2lp (
    dsp_zpk *zpk,
    double   Wc
)
{
    dsp_zpk *zpk_lp;
    CALLOC_RET(zpk_lp, 1, dsp_zpk, NULL);

    zpk_lp->n_zero = zpk->n_zero;
    zpk_lp->n_pole = zpk->n_pole;
    
    if( zpk_lp->n_zero > 0 ){
        CALLOC_GOTO(zpk_lp->zero, zpk_lp->n_zero, doubleC, FALSE1);

        for(uint16_t i = 0; i < zpk_lp->n_zero; i++){
            zpk_lp->zero[i] = Wc*zpk->zero[i];
        }
    } else {
        zpk_lp->zero = NULL;
    }

    if( zpk_lp->n_pole > 0 ){
        CALLOC_GOTO(zpk_lp->pole, zpk_lp->n_pole, doubleC, FALSE2);

        for(uint16_t i = 0; i < zpk_lp->n_pole; i++){
            zpk_lp->pole[i] = Wc*zpk->pole[i];
        }
    } else {
        zpk_lp->pole = NULL;
    }

    zpk_lp->gain = zpk->gain*pow(Wc, (zpk->n_pole - zpk->n_zero));
    return zpk_lp;

FALSE2:
    if( zpk_lp->pole )  free(zpk_lp->pole);
FALSE1:
    if( zpk_lp )    free(zpk_lp);
    return NULL;
}

dsp_zpk *zpklp2hp (
    dsp_zpk *zpk,
    double   Wc
)
{
    uint16_t degree;
    degree = zpk->n_pole - zpk->n_zero;
    if( degree < 0 ){
        fprintf(stderr, "Error: Not stable filter\n");
        return NULL;
    }

    dsp_zpk *zpk_hp;
    CALLOC_RET(zpk_hp, 1, dsp_zpk, NULL);

    zpk_hp->n_zero = zpk->n_pole;
    zpk_hp->n_pole = zpk->n_pole;

    if( zpk->n_pole == 0 ){
        zpk_hp->zero = NULL;
        zpk_hp->pole = NULL;
        return zpk_hp;
    }

    uint16_t i;
    doubleC  gain = 1.0;
    CALLOC_GOTO(zpk_hp->zero, zpk_hp->n_zero, doubleC, FALSE1);
    for(i = 0; i <    zpk->n_zero; i++){
        zpk_hp->zero[i] = Wc/zpk->zero[i];
        gain *= (-zpk->zero[i]);
    }
    for(     ; i < zpk_hp->n_zero; i++){
        zpk_hp->zero[i] = 0.0;
    }

    CALLOC_GOTO(zpk_hp->pole, zpk_hp->n_pole, doubleC, FALSE2);
    for(i = 0; i < zpk_hp->n_pole; i++){
        zpk_hp->pole[i] = Wc/zpk->pole[i];
        gain /= (-zpk->pole[i]);
    }

    zpk_hp->gain = zpk->gain*creal(gain);
    return zpk_hp;

FALSE2:
    if( zpk_hp->pole )  free(zpk_hp->pole);
FALSE1:
    if( zpk_hp )    free(zpk_hp);
    return NULL;
}

dsp_zpk *zpklp2bp (
    dsp_zpk *zpk,
    double   Wc[]
)
{
    double bw = (Wc[1] - Wc[0]);
    double Wo = sqrt( Wc[0]*Wc[1] );

    uint16_t degree;
    degree = zpk->n_pole - zpk->n_zero;
    if( degree < 0 ){
        fprintf(stderr, "Error: Not stable filter\n");
        return NULL;
    }

    dsp_zpk *zpk_bp;
    CALLOC_RET(zpk_bp, 1, dsp_zpk, NULL);;

    doubleC *z_lp, *p_lp;
    CALLOC_GOTO(z_lp, zpk->n_zero, doubleC, FALSE1);
    CALLOC_GOTO(p_lp, zpk->n_pole, doubleC, FALSE2);

    for(uint16_t i = 0; i < zpk->n_zero; i++){
        z_lp[i] = zpk->zero[i] * bw/2.0;
    }
    for(uint16_t i = 0; i < zpk->n_pole; i++){
        p_lp[i] = zpk->pole[i] * bw/2.0;
    }

    zpk_bp->n_zero = zpk->n_zero*2 + degree;
    zpk_bp->n_pole = zpk->n_pole*2;

    if( zpk_bp->n_pole > 0 ){
        CALLOC_GOTO(zpk_bp->pole, zpk_bp->n_pole, doubleC, FALSE3);
        for(uint16_t i = 0; i < zpk->n_pole; i++){
            uint16_t j = 2*i;
            uint16_t k = 2*i + 1;
            zpk_bp->pole[j] = p_lp[i] + csqrt(p_lp[i]*p_lp[i] - Wo*Wo);
            zpk_bp->pole[k] = p_lp[i] - csqrt(p_lp[i]*p_lp[i] - Wo*Wo);
#ifdef DEBUG
            fprintf(stderr, "pole_bp[%u] = %e + I*%e\n", j, creal(zpk_bp->pole[j]), cimag(zpk_bp->pole[j]));
            fprintf(stderr, "pole_bp[%u] = %e + I*%e\n", k, creal(zpk_bp->pole[k]), cimag(zpk_bp->pole[k]));
#endif
        }
    } else {
        zpk_bp->pole = NULL;
    }

    if( zpk_bp->n_zero > 0 ){
        CALLOC_GOTO(zpk_bp->zero, zpk_bp->n_zero, doubleC, FALSE4);
        for(uint16_t i = 0; i < zpk->n_zero; i++){
            uint16_t j = 2*i;
            uint16_t k = 2*i + 1;
            zpk_bp->zero[j] = z_lp[i] + csqrt(z_lp[i]*z_lp[i] - Wo*Wo);
            zpk_bp->zero[k] = z_lp[i] - csqrt(z_lp[i]*z_lp[i] - Wo*Wo);
#ifdef DEBUG
            fprintf(stderr, "zero_bp[%u] = %e + I*%e\n", j, creal(zpk_bp->zero[j]), cimag(zpk_bp->zero[j]));
            fprintf(stderr, "zero_bp[%u] = %e + I*%e\n", k, creal(zpk_bp->zero[k]), cimag(zpk_bp->zero[k]));
#endif
        }

        for(uint16_t i = 0; i < degree; i++){
            uint16_t j = i + zpk->n_zero*2;
            zpk_bp->zero[j] = 0;
#ifdef DEBUG
            fprintf(stderr, "zero_bp[%u] = %e + I*%e\n", j, creal(zpk_bp->zero[j]), cimag(zpk_bp->zero[j]));
#endif
        }
    } else {
        zpk_bp->zero = NULL;
    }

    zpk_bp->gain = zpk->gain*pow(bw, degree);

    free(z_lp);
    free(p_lp);

    return zpk_bp;

FALSE4:
    if( zpk_bp->pole )  free(zpk_bp->pole);
FALSE3:
    free(p_lp);
FALSE2:
    free(z_lp);
FALSE1:
    if( zpk_bp )    free(zpk_bp);
    return NULL;
}

dsp_zpk *zpklp2bs (
    dsp_zpk *zpk,
    double   Wc[]
)
{
    double bw = (Wc[1] - Wc[0]);
    double Wo = sqrt( Wc[0]*Wc[1] );

    uint16_t degree;
    degree = zpk->n_pole - zpk->n_zero;
    if( degree < 0 ){
        fprintf(stderr, "Error: Not stable filter\n");
        return NULL;
    }

    dsp_zpk *zpk_bs;
    CALLOC_RET(zpk_bs, 1, dsp_zpk, NULL);;

    doubleC *z_hp, *p_hp;
    CALLOC_GOTO(z_hp, zpk->n_zero, doubleC, FALSE1);
    CALLOC_GOTO(p_hp, zpk->n_pole, doubleC, FALSE2);

    for(uint16_t i = 0; i < zpk->n_zero; i++){
        z_hp[i] = (bw/2.0) / zpk->zero[i];
    }
    for(uint16_t i = 0; i < zpk->n_pole; i++){
        p_hp[i] = (bw/2.0) / zpk->pole[i];
    }

    zpk_bs->n_zero = zpk->n_zero*2 + degree*2;
    zpk_bs->n_pole = zpk->n_pole*2;

    doubleC  gain = 1.0;
    if( zpk_bs->n_pole > 0 ){
        CALLOC_GOTO(zpk_bs->pole, zpk_bs->n_pole, doubleC, FALSE3);
        for(uint16_t i = 0; i < zpk->n_pole; i++){
            uint16_t j = 2*i;
            uint16_t k = 2*i + 1;
            zpk_bs->pole[j] = p_hp[i] + csqrt(p_hp[i]*p_hp[i] - Wo*Wo);
            zpk_bs->pole[k] = p_hp[i] - csqrt(p_hp[i]*p_hp[i] - Wo*Wo);
            gain *= zpk_bs->pole[j]*zpk_bs->pole[k];
#ifdef DEBUG
            fprintf(stderr, "pole_bs[%u] = %e + I*%e\n", j, creal(zpk_bs->pole[j]), cimag(zpk_bs->pole[j]));
            fprintf(stderr, "pole_bs[%u] = %e + I*%e\n", k, creal(zpk_bs->pole[k]), cimag(zpk_bs->pole[k]));
#endif
        }
    } else {
        zpk_bs->pole = NULL;
        gain = 0.0;
    }

    if( zpk_bs->n_zero > 0 ){
        CALLOC_GOTO(zpk_bs->zero, zpk_bs->n_zero, doubleC, FALSE4);
        for(uint16_t i = 0; i < zpk->n_zero; i++){
            uint16_t j = 2*i;
            uint16_t k = 2*i + 1;
            zpk_bs->zero[j] = z_hp[i] + csqrt(z_hp[i]*z_hp[i] - Wo*Wo);
            zpk_bs->zero[k] = z_hp[i] - csqrt(z_hp[i]*z_hp[i] - Wo*Wo);
            gain /= zpk_bs->zero[j]*zpk_bs->zero[k];
#ifdef DEBUG
            fprintf(stderr, "zero_bs[%u] = %e + I*%e\n", j, creal(zpk_bs->zero[j]), cimag(zpk_bs->zero[j]));
            fprintf(stderr, "zero_bs[%u] = %e + I*%e\n", k, creal(zpk_bs->zero[k]), cimag(zpk_bs->zero[k]));
#endif
        }

        for(uint16_t i = 0; i < degree; i++){
            uint16_t j = zpk->n_zero*2 + i*2;
            uint16_t k = zpk->n_zero*2 + i*2 + 1;
            zpk_bs->zero[j] =  I*Wo;
            zpk_bs->zero[k] = -I*Wo;
            gain /= zpk_bs->zero[j]*zpk_bs->zero[k];
#ifdef DEBUG
            fprintf(stderr, "zero_bs[%u] = %e + I*%e\n", j, creal(zpk_bs->zero[j]), cimag(zpk_bs->zero[j]));
            fprintf(stderr, "zero_bs[%u] = %e + I*%e\n", k, creal(zpk_bs->zero[k]), cimag(zpk_bs->zero[k]));
#endif
        }
    } else {
        zpk_bs->zero = NULL;
    }

    zpk_bs->gain = zpk->gain*creal(gain);

    free(z_hp);
    free(p_hp);

    return zpk_bs;

FALSE4:
    if( zpk_bs->pole )  free(zpk_bs->pole);
FALSE3:
    free(p_hp);
FALSE2:
    free(z_hp);
FALSE1:
    if( zpk_bs )    free(zpk_bs);
    return NULL;
}


dsp_zpk *zpk_bilinear (
    dsp_zpk *zpk
)
{
    if( zpk->n_zero > zpk->n_pole ){
        fprintf(stderr, "must be n_zero <= n_pole\n");
        return NULL;
    }

    dsp_zpk *zpk_z;
    CALLOC_GOTO(zpk_z, 1, dsp_zpk, FALSE1);

    /* size of n_zero is equal to n_pole by fixing -1 */
    zpk_z->n_zero = zpk->n_pole;
    zpk_z->n_pole = zpk->n_pole;

    CALLOC_GOTO(zpk_z->zero, zpk_z->n_zero, doubleC, FALSE2);
    CALLOC_GOTO(zpk_z->pole, zpk_z->n_zero, doubleC, FALSE3);

    uint16_t i;
    double  fs2 = 4.0;
    for(i = 0; i < zpk->n_zero; i++){
        zpk_z->zero[i] = (fs2 + zpk->zero[i])/(fs2 - zpk->zero[i]);
#ifdef DEBUG
        fprintf(stderr, "zero_z[%u] = %e + I*%e\n", i, creal(zpk_z->zero[i]), cimag(zpk_z->zero[i]));
#endif
    }
    for(     ; i < zpk_z->n_zero; i++){
        zpk_z->zero[i] = -1;
#ifdef DEBUG
        fprintf(stderr, "zero_z[%u] = %e + I*%e\n", i, creal(zpk_z->zero[i]), cimag(zpk_z->zero[i]));
#endif
    }

    for(i = 0; i < zpk->n_pole; i++){
        zpk_z->pole[i] = (fs2 + zpk->pole[i])/(fs2 - zpk->pole[i]);
#ifdef DEBUG
        fprintf(stderr, "pole_z[%u] = %e + I*%e\n", i, creal(zpk_z->pole[i]), cimag(zpk_z->pole[i]));
#endif
    }

    doubleC  gain = 1.0;
    for(i = 0; i < zpk->n_zero; i++){
        gain *= (fs2 - zpk->zero[i])/(fs2 - zpk->pole[i]);
    }
    for(     ; i < zpk->n_pole; i++){
        gain /= (fs2 - zpk->pole[i]);
    }
    zpk_z->gain = zpk->gain*creal(gain);
#ifdef DEBUG
    fprintf(stderr, "gain_z = %e\n", zpk_z->gain);
#endif

    return zpk_z;

FALSE3:
    if( zpk_z->zero)    free(zpk_z->zero);
FALSE2:
    if( zpk_z ) free(zpk_z);
FALSE1:
    return NULL;
}


dsp_tf *zpk2tf (
    dsp_zpk *zpk
)
{
    doubleC *a = NULL;
    doubleC *b = NULL;

    a = poly(zpk->pole, zpk->n_pole);
    b = poly(zpk->zero, zpk->n_zero);
    if(( a == NULL )||( b == NULL ))    goto FALSE1;

    dsp_tf  *tf = NULL;
    CALLOC_GOTO(tf, 1, dsp_tf, FALSE2);

    tf->n_a = zpk->n_pole + 1;
    tf->n_b = zpk->n_zero + 1;
    CALLOC_GOTO(tf->a, tf->n_a, double, FALSE2);
    CALLOC_GOTO(tf->b, tf->n_b, double, FALSE3);

    printf("gain = %e\n", zpk->gain);
    for(uint16_t i = 0; i < tf->n_a; i++){
        tf->a[i] = creal(a[i]);
//      fprintf(stderr, "a[%u] = %.10e + I*%.10e\n", i, creal(a[i]), cimag(a[i]));
#ifdef DEBUG
        fprintf(stderr, "a[%u] = %.10e\n", i, tf->a[i]);
#endif
    }
    for(uint16_t i = 0; i < tf->n_b; i++){
        tf->b[i] = zpk->gain*creal(b[i]);
//      fprintf(stderr, "b[%u] = %.10e + I*%.10e\n", i, creal(b[i]), cimag(b[i]));
#ifdef DEBUG
        fprintf(stderr, "b[%u] = %.10e\n", i, tf->b[i]);
#endif
    }

    free(a);
    free(b);
    return tf;

FALSE3:
    free(tf->a);
FALSE2:
    free(tf);
FALSE1:
    if(b)   free(b);
    if(a)   free(a);
    return NULL;
}

dsp_tf *bilinear (
    dsp_tf *tf_s
)
{
    dsp_tf *tf;
    CALLOC_RET(tf, 1, dsp_tf, NULL);

    tf->n_a = tf_s->n_a;
    tf->n_b = tf_s->n_a;

    CALLOC_GOTO(tf->a, tf->n_a, double, FALSE1);
    CALLOC_GOTO(tf->b, tf->n_b, double, FALSE2);

    for(uint16_t k = 0; k < tf->n_a; k++){
        double a = 1;

        for(uint16_t j = 0; j < (tf->n_a - k); j++){
            a += nCr(tf->n_a, k)*tf_s->a[j];
        }
        if( k&1 )   a *= (-1);
        tf->a[k] = a;
#ifdef DEBUG
        fprintf(stderr, "a_z[%u] = %e\n", k, tf->a[k]);
#endif
    }

    for(uint16_t n = 0; n < tf->n_b; n++){
        double b = 1;

        for(uint16_t l = 0; l <= n; l++){
            double bb = 0;

            for(uint16_t i = 0; i < (tf_s->n_b - l); i++){
                bb += nCr((tf_s->n_b - i), l)*tf_s->b[i];
            }
            bb *= nCr((tf_s->n_a - tf_s->n_b), (n - l));
            if( l&1 )   bb *= (-1);
            b += bb;
        }
        b *= pow(2, (tf_s->n_b - tf_s->n_a));
        tf->b[n] = b;
#ifdef DEBUG
        fprintf(stderr, "b_z[%u] = %e\n", n, tf->b[n]);
#endif
    }

FALSE2:
    free(tf->a);
FALSE1:
    free(tf);
    return NULL;
}

/**
 * Get the second-order-sections representation of the digital filter with the specified zero-pole-gain
 * zpk: zero-pole-gain of the filter
 */
dsp_sos *zpk2sos (
    dsp_zpk  *zpk

)
{
    if( zpk->n_zero > zpk->n_pole ){
        fprintf(stderr, "Not stable filter.\n");
        return NULL;
    }

    dsp_sos *sos;
    CALLOC_RET(sos, 1, dsp_sos, NULL);

    uint16_t L = (zpk->n_pole + 1)/2;
    sos->L = L;
    CALLOC_GOTO(sos->a, L, double*, FALSE1);
    for(uint16_t i = 0; i < L; i++){
        CALLOC_GOTO(sos->a[i], 6, double, FALSE2);
    }
    
    doubleC *z, *p;
    doubleC *z_neg, *p_neg;
    doubleC *z_sort, *p_sort;
    CALLOC_GOTO(z,      zpk->n_zero,   doubleC, FALSE2);
    CALLOC_GOTO(p,      zpk->n_pole,   doubleC, FALSE3);
    CALLOC_GOTO(z_neg,  zpk->n_zero/2, doubleC, FALSE4);
    CALLOC_GOTO(p_neg,  zpk->n_pole/2, doubleC, FALSE5);
    CALLOC_GOTO(z_sort, zpk->n_pole,   doubleC, FALSE6);
    CALLOC_GOTO(p_sort, zpk->n_pole,   doubleC, FALSE7);
    
    for(uint16_t i = 0; i < zpk->n_zero; i++){
        z[i] = zpk->zero[i];
    }
    for(uint16_t i = 0; i < zpk->n_pole; i++){
        p[i] = zpk->pole[i];
    }

    /* When the order is odd, there is one real pole/zero.
       It should be ignored in making pair.
       Here it moves to the last */
    if( zpk->n_zero%2 == 1 ){
        double tmp = z[zpk->n_zero/2];
        for(uint16_t i = zpk->n_zero/2; i < (zpk->n_zero - 1); i++){
            z[i] = z[i+1];
        }
        z[zpk->n_zero - 1] = tmp;
    }
    if( zpk->n_pole%2 == 1 ){
        double tmp = p[zpk->n_pole/2];
        for(uint16_t i = zpk->n_pole/2; i < (zpk->n_pole - 1); i++){
            p[i] = p[i+1];
        }
        p[zpk->n_pole - 1] = tmp;
    }

    /* When the order is odd, ignore the last element */
    cplxpair(z, zpk->n_zero - zpk->n_zero%2);
    cplxpair(p, zpk->n_pole - zpk->n_pole%2);

    //cplxpair(z, zpk->n_zero);
    //cplxpair(p, zpk->n_pole);

#if defined( DEBUG )
    {
        for(uint16_t i = 0; i < zpk->n_zero; i++){
            fprintf(stderr, "z[%u] = %e + I*%e\n", i, creal(z[i]), cimag(z[i]));
        }
        for(uint16_t i = 0; i < zpk->n_pole; i++){
        fprintf(stderr, "p[%u] = %e + I*%e\n", i, creal(p[i]), cimag(p[i]));
        }
    }
#endif
    
    for(uint16_t i = 0; i < zpk->n_zero/2; i++){
        z_neg[i] = z[2*i];
    }
    for(uint16_t i = 0; i < zpk->n_pole/2; i++){
        p_neg[i] = p[2*i];
    }
    
    uint16_t  i_min = 0;
    double    r_min = 1e5;
    uint8_t  *picked;

    CALLOC_GOTO(picked, zpk->n_zero/2, uint8_t, FALSE8);
    
    for(uint16_t i = 0; i < zpk->n_zero/2; i++){
        for(uint16_t j = 0; j < zpk->n_zero/2; j++){
            if( picked[j] )    continue;

            double r = cabs(p_neg[i] - z_neg[j]);
            if( r < r_min ){
                r_min = r;
                i_min = j;
            }
        }
        
        p_sort[2*i]   = p[2*i];
        p_sort[2*i+1] = p[2*i+1];
        z_sort[2*i]   = z[2*i_min];
        z_sort[2*i+1] = z[2*i_min+1];
        picked[i_min] = 1;
        r_min = 1e5;
    }

    for(uint16_t i = zpk->n_zero; i < zpk->n_pole; i++){
      z_sort[i] = 0.0;
      p_sort[i] = p[i];
    }
    
    //when N is odd, say N = 5, z[5] and p[5] are undefined!!
    if( zpk->n_zero%2 == 1 ){
        z_sort[zpk->n_zero - 1] = z[zpk->n_zero - 1];
    }
    if( zpk->n_pole%2 == 1 ){
        p_sort[zpk->n_pole - 1] = p[zpk->n_pole - 1];
    }

    for(uint16_t i = 0; i < L; i++){
        uint16_t j = 2*i;
        uint16_t k = 2*i + 1;

        sos->a[i][0] = 1.0;
        sos->a[i][1] = -(z_sort[j] + z_sort[k]);
        sos->a[i][2] = z_sort[j]*z_sort[k];
        sos->a[i][3] = 1.0;
        sos->a[i][4] = -(p_sort[j] + p_sort[k]);
        sos->a[i][5] = p_sort[j]*p_sort[k];

#ifdef DEBUG
        for(uint8_t m = 0; m < 6; m++){
            fprintf(stderr, " %.5e", sos->a[i][m]);
        }
        fprintf(stderr, "\n");
#endif
    }
    sos->G = zpk->gain;
    
    free(p);
    free(z);
    free(z_neg);
    free(p_neg);
    free(z_sort);
    free(p_sort);
    free(picked);
    return sos;


    free(picked);
FALSE8:
    free(p_sort);
FALSE7:
    free(z_sort);
FALSE6:
    free(p_neg);
FALSE5:
    free(z_neg);
FALSE4:
    free(p);
FALSE3:
    free(z);
FALSE2:
    for(uint16_t i = 0; i < L; i++){
        if( sos->a[i] )    free(sos->a[i]);
    }
    free(sos->a);
FALSE1:
    free(sos);
    return NULL;
}

dsp_sos *iir_design_sos (
    double   wp[],
    double   ws[],
    double   Gp,
    double   As,
    iir_type type
)
{
    // uint8_t     n_f;
    filter_type ftype;

    /* lowpass or highpass */
    if(( wp[1] < 0 )||( ws[1] < 0 )){
        /* lowpass */
        if( wp[0] <= ws[0] ){
            // n_f = 1;
            ftype = LOWPASS;
        }
        /* highpass */
        else{
            // n_f = 1;
            ftype = HIGHPASS;
        }
    }
    else{
        if( wp[0] >= wp[1] ){
            fprintf(stderr, "to bandpass/bandstop: wp[0] < wp[1]\n");
            return NULL;
        }
        if( ws[0] >= ws[1] ){
            fprintf(stderr, "to bandpass/bandstop: ws[0] < ws[1]\n");
            return NULL;
        }

        /* bandpass */
        if( wp[0] >= ws[0] ){
            if( wp[1] > ws[1] ){
                fprintf(stderr, "to bandpass: ws[0] < wp[0] < wp[1] < ws[1]\n");
                return NULL;
            }

            // n_f = 2;
            ftype = BANDPASS;
        }
        /* bandstop */
        else{
            if( wp[1] < ws[1] ){
                fprintf(stderr, "to bandstop: wp[0] < ws[0] < ws[1] < wp[1]\n");
                return NULL;
            }

            // n_f = 2;
            ftype = BANDSTOP;
        }
    }

    dsp_butter_arg *butarg;
    butarg = butter_order(wp, ws, Gp, As, ftype);

    butarg->order += 1;
    fprintf(stderr, "order = %u\n", butarg->order);

#ifdef DEBUG
    fprintf(stderr, "ws[0] = %f\n", ws[0]);
    fprintf(stderr, "wc[0] = %f\n", butarg->wc[0]);
    fprintf(stderr, "wp[0] = %f\n", wp[0]);
    fprintf(stderr, "wp[1] = %f\n", wp[1]);
    fprintf(stderr, "wc[1] = %f\n", butarg->wc[1]);
    fprintf(stderr, "ws[1] = %f\n", ws[1]);
#endif

    dsp_zpk *zpk, *zpk2, *zpk_dig;
    zpk = butter_zpk(butarg->order);
    
    /* Pre-warping */
    double Wc[2] = {0};
    double fs = 2.0;
    Wc[0] = 2.0*fs*tan(M_PI*butarg->wc[0]/fs);
    if(( ftype == BANDPASS )||( ftype == BANDSTOP )){
        Wc[1] = 2.0*fs*tan(M_PI*butarg->wc[1]/fs);
    }

#ifdef DEBUG
    fprintf(stderr, "Wc[0] = %f\n", Wc[0]);
    fprintf(stderr, "Wc[1] = %f\n", Wc[1]);
#endif

    switch( ftype ){
        case LOWPASS:
            zpk2 = zpklp2lp(zpk, Wc[0]);
            break;
        case HIGHPASS:
            zpk2 = zpklp2hp(zpk, Wc[0]);
            break;
        case BANDPASS:
            zpk2 = zpklp2bp(zpk, Wc);
            break;
        case BANDSTOP:
            zpk2 = zpklp2bs(zpk, Wc);
            break;
        default:
            return NULL;
    }

/*
    dsp_tf  *tf_a;
    tf_a = zpk2tf(zpk2);

    dsp_tf  *tf_z;
    tf_z = bilinear(tf_a);
    exit(0);
    */
//  fprintf(stderr, "gain: %e\n", zpk2->gain);
/*
    {
        freq_res_s  *fres;
        //fres = freq_responce_s(zpk2, NULL, 2048);
        fres = freq_responce_s(zpk, NULL, 2048);
        print_fres_s(fres);
    }
    exit(0);
*/
    zpk_dig = zpk_bilinear(zpk2);
//  fprintf(stderr, "gain_dig: %e\n", zpk_dig->gain);
/*
    {
        freq_res_z  *fres;
        fres = freq_responce_z(zpk_dig, NULL, 8192);
        print_fres_z(fres);
    }
    exit(0);
*/
/*
    dsp_tf  *tf;
    tf = zpk2tf(zpk_dig);

    exit(0);
    {
        freq_res_z  *fres;
        fres = freq_responce_z(NULL, tf, 1024);
        print_fres_z(fres);
    }
    exit(0);
*/

    dsp_sos  *sos;
    sos = zpk2sos(zpk_dig);

    free(zpk);
    free(zpk2);
    free(zpk_dig);

    return sos;
}

void free_dsp_tf (
    dsp_tf  *tf
)
{
    if( tf->a ) free(tf->a);
    if( tf->b ) free(tf->b);
    free(tf);
}

void free_dsp_sos (
    dsp_sos  *sos
)
{
    for(uint16_t i = 0; i < sos->L; i++){
        free(sos->a[i]);
    }
    free(sos->a);
    free(sos);
}
