#include <stdio.h>
#include <stdlib.h>

void polint(const double *xa, const double *ya, double x,
            double *y, double *dy, int ordn);

void polin3(const double *x1a, const double *x2a, const double *x3a,
            const double *ya, double x1, double x2, double x3,
            double *y, double *dy, int ordn)
{
    #define YA3(i,j,k) ya[(i) + ordn*((j) + ordn*(k))]

    int j, k;
    double dy_temp;

    double *yatmp = (double*)malloc((size_t)ordn * (size_t)ordn * sizeof(double));
    double *ymtmp = (double*)malloc((size_t)ordn * sizeof(double));
    if (!yatmp || !ymtmp) {
        fprintf(stderr, "polin3: malloc failed\n");
        exit(1);
    }
    #define YAT(j,k) yatmp[(j) + ordn*(k)]

    for (k = 0; k < ordn; k++) {
        for (j = 0; j < ordn; j++) {
            polint(x1a, &YA3(0, j, k), x1, &YAT(j, k), &dy_temp, ordn);
        }
    }

    for (k = 0; k < ordn; k++) {
        polint(x2a, &YAT(0, k), x2, &ymtmp[k], &dy_temp, ordn);
    }

    polint(x3a, ymtmp, x3, y, dy, ordn);

    #undef YAT
    free(yatmp);
    free(ymtmp);
    #undef YA3
}
