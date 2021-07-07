#include <math.h>
#include <stdio.h>
#include "dsplib.h"

int main(void)
{
    double fs = 4096;
    int k;
    double t, x;

    for(k = 0; k < 1024; k++){
        t = k/fs;
        x = sin(80*M_PIl*t) + 0.3*sin(300*M_PIl*t) + 0.7*sin(700*M_PI*t);
        printf("%.10e %.10e\n", t, x);
    }

    return 0;
}
