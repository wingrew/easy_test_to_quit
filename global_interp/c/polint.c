// polint_polin23.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
  Fortran layout (column-major):
    ya(i,j)   -> ya[(i-1) + ordn*(j-1)]
    ya(i,j,k) -> ya[(i-1) + ordn*((j-1) + ordn*(k-1))]
  In C (0-based):
    YA2(i,j)  -> ya[(i) + ordn*(j)]
    YA3(i,j,k)-> ya[(i) + ordn*((j) + ordn*(k))]
*/

void polint(const double *xa, const double *ya, double x,
                   double *y, double *dy, int ordn)
{
    int i, m, ns, n_m;
    double dif, dift, hp, h, den_val;

    double *c  = (double*)malloc((size_t)ordn * sizeof(double));
    double *d  = (double*)malloc((size_t)ordn * sizeof(double));
    double *ho = (double*)malloc((size_t)ordn * sizeof(double));
    if (!c || !d || !ho) {
        fprintf(stderr, "polint: malloc failed\n");
        exit(1);
    }

    for (i = 0; i < ordn; i++) {
        c[i]  = ya[i];
        d[i]  = ya[i];
        ho[i] = xa[i] - x;
    }

    ns  = 0;                      // Fortran ns=1 -> C ns=0
    dif = fabs(x - xa[0]);

    for (i = 1; i < ordn; i++) {
        dift = fabs(x - xa[i]);
        if (dift < dif) {
            ns  = i;
            dif = dift;
        }
    }

    *y  = ya[ns];
    ns -= 1;                      // Fortran ns=ns-1

    for (m = 1; m <= ordn - 1; m++) {
        n_m = ordn - m;           // number of active points this round
        for (i = 0; i < n_m; i++) {
            hp      = ho[i];
            h       = ho[i + m];
            den_val = hp - h;

            if (den_val == 0.0) {
                fprintf(stderr, "failure in polint for point %g\n", x);
                fprintf(stderr, "with input points xa: ");
                for (int t = 0; t < ordn; t++) fprintf(stderr, "%g ", xa[t]);
                fprintf(stderr, "\n");
                exit(1);
            }

            den_val = (c[i + 1] - d[i]) / den_val;
            d[i]    = h  * den_val;
            c[i]    = hp * den_val;
        }

        // Fortran: if (2*ns < n_m) then dy=c(ns+1) else dy=d(ns); ns=ns-1
        // Here ns is C-indexed and can be -1; logic still matches.
        if (2 * ns < n_m) {
            *dy = c[ns + 1];
        } else {
            *dy = d[ns];
            ns -= 1;
        }
        *y += *dy;
    }

    free(c);
    free(d);
    free(ho);
}



