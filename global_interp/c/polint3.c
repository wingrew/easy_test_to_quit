void polin3(const double *x1a, const double *x2a, const double *x3a,
                   const double *ya, double x1, double x2, double x3,
                   double *y, double *dy, int ordn)
{
    // ya is ordn x ordn x ordn in Fortran layout (column-major)
    #define YA3(i,j,k) ya[(i) + ordn*((j) + ordn*(k))]  // i,j,k: 0..ordn-1

    int j, k;
    double dy_temp;

    // yatmp(j,k) in Fortran code is ordn x ordn, treat column-major:
    // yatmp(j,k) -> yatmp[j + ordn*k]
    double *yatmp = (double*)malloc((size_t)ordn * (size_t)ordn * sizeof(double));
    double *ymtmp = (double*)malloc((size_t)ordn * sizeof(double));
    if (!yatmp || !ymtmp) {
        fprintf(stderr, "polin3: malloc failed\n");
        exit(1);
    }
    #define YAT(j,k) yatmp[(j) + ordn*(k)]

    for (k = 0; k < ordn; k++) {
        for (j = 0; j < ordn; j++) {
            // call polint(x1a, ya(:,j,k), x1, yatmp(j,k), dy_temp)
            // ya(:,j,k) contiguous: base is &YA3(0,j,k)
            polint(x1a, &YA3(0, j, k), x1, &YAT(j, k), &dy_temp, ordn);
        }
    }

    for (k = 0; k < ordn; k++) {
        // call polint(x2a, yatmp(:,k), x2, ymtmp(k), dy_temp)
        polint(x2a, &YAT(0, k), x2, &ymtmp[k], &dy_temp, ordn);
    }

    polint(x3a, ymtmp, x3, y, dy, ordn);

    #undef YAT
    free(yatmp);
    free(ymtmp);
    #undef YA3
}