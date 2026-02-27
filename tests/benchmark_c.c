#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void global_interp(const int ex[3],
                   const double *X, const double *Y, const double *Z,
                   const double *f,
                   double *f_int,
                   double x1, double y1, double z1,
                   int ORDN,
                   const double SoA[3],
                   int symmetry);

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s <num_queries> <output_file>\n", argv[0]);
        return 2;
    }

    const int nquery = atoi(argv[1]);
    const char *out_path = argv[2];

    const int ex[3] = {24, 20, 18};
    const int ordn = 4;
    const int symmetry = 2;
    const double soa[3] = {1.0, -1.0, 1.0};

    const double dx = 0.1;
    const double dy = 0.2;
    const double dz = 0.15;

    double *X = (double *)malloc((size_t)ex[0] * sizeof(double));
    double *Y = (double *)malloc((size_t)ex[1] * sizeof(double));
    double *Z = (double *)malloc((size_t)ex[2] * sizeof(double));
    double *f = (double *)malloc((size_t)ex[0] * (size_t)ex[1] * (size_t)ex[2] * sizeof(double));

    if (!X || !Y || !Z || !f) {
        fprintf(stderr, "alloc failed\n");
        return 2;
    }

    for (int i = 0; i < ex[0]; ++i) X[i] = i * dx;
    for (int j = 0; j < ex[1]; ++j) Y[j] = j * dy;
    for (int k = 0; k < ex[2]; ++k) Z[k] = k * dz;

    for (int k = 1; k <= ex[2]; ++k) {
        for (int j = 1; j <= ex[1]; ++j) {
            for (int i = 1; i <= ex[0]; ++i) {
                const double x = X[i - 1];
                const double y = Y[j - 1];
                const double z = Z[k - 1];
                const size_t idx = (size_t)(i - 1) + (size_t)ex[0] * ((size_t)(j - 1) + (size_t)ex[1] * (size_t)(k - 1));
                f[idx] = sin(x) + 0.5 * cos(2.0 * y) + z * z + 0.1 * x * y * z;
            }
        }
    }

    FILE *fp = fopen(out_path, "w");
    if (!fp) {
        perror("fopen");
        return 2;
    }

    uint32_t state = 123456789u;
    const double x_max = X[ex[0] - 1];
    const double y_max = Y[ex[1] - 1];
    const double z_max = Z[ex[2] - 1];

    double t0 = now_sec();
    double checksum = 0.0;
    for (int q = 0; q < nquery; ++q) {
        state = 1664525u * state + 1013904223u;
        const double rx = (double)state / (double)UINT32_MAX;
        state = 1664525u * state + 1013904223u;
        const double ry = (double)state / (double)UINT32_MAX;
        state = 1664525u * state + 1013904223u;
        const double rz = (double)state / (double)UINT32_MAX;

        const double xq = rx * x_max;
        const double yq = ry * y_max;
        const double zq = rz * z_max;

        double out;
        global_interp(ex, X, Y, Z, f, &out, xq, yq, zq, ordn, soa, symmetry);
        checksum += out;
        fprintf(fp, "%.17g\n", out);
    }
    double elapsed = now_sec() - t0;

    fclose(fp);

    fprintf(stderr, "elapsed_sec=%.9f\n", elapsed);
    fprintf(stderr, "checksum=%.17g\n", checksum);

    free(X);
    free(Y);
    free(Z);
    free(f);
    return 0;
}
