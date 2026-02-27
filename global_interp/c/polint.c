#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void polint(const double *xa, const double *ya, double x,
            double *y, double *dy, int ordn)
{
    int i, m, n_m;
    int ns_f;
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

    ns_f = 1;
    dif = fabs(x - xa[0]);

    for (i = 1; i < ordn; i++) {
        dift = fabs(x - xa[i]);
        if (dift < dif) {
            ns_f = i + 1;
            dif = dift;
        }
    }

    *y = ya[ns_f - 1];
    ns_f -= 1;

    for (m = 1; m <= ordn - 1; m++) {
        n_m = ordn - m;
        for (i = 0; i < n_m; i++) {
            hp = ho[i];
            h = ho[i + m];
            den_val = hp - h;

            if (den_val == 0.0) {
                fprintf(stderr, "failure in polint for point %g\n", x);
                fprintf(stderr, "with input points xa: ");
                for (int t = 0; t < ordn; t++) fprintf(stderr, "%g ", xa[t]);
                fprintf(stderr, "\n");
                exit(1);
            }

            den_val = (c[i + 1] - d[i]) / den_val;
            d[i] = h * den_val;
            c[i] = hp * den_val;
        }

        if (2 * ns_f < n_m) {
            *dy = c[ns_f];
        } else {
            *dy = d[ns_f - 1];
            ns_f -= 1;
        }
        *y += *dy;
    }

    free(c);
    free(d);
    free(ho);
}
