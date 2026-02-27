#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 你已有的 polin3（由前面 Fortran->C 翻译得到） */
void polin3(const double *x1a, const double *x2a, const double *x3a,
            const double *ya, double x1, double x2, double x3,
            double *y, double *dy, int ordn);

/*
  你需要提供 decide3d 的实现（这里仅声明）。
  Fortran: decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry)
  - ex: [3]
  - f: 三维场（列主序）
  - cxB/cxT: 3 维窗口起止（Fortran 1-based，且可能 <=0）
  - SoA: [3]
  - ya: 输出 ORDN^3 的采样块（列主序）
  - return: 0 表示正常；非 0 表示错误（对应 Fortran logical = .true.）
*/
int decide3d(const int ex[3],
             const double *f_in,
             const double *f_in2,   /* Fortran 里传了 f,f；按原样保留 */
             const int cxB[3],
             const int cxT[3],
             const double SoA[3],
             double *ya,
             int ordn,
             int symmetry);

/* 把 Fortran 1-based 下标 idxF (可为负/0) 映射到 C 的 X[idx] 访问（只用于 X(2-cxB) 这种表达式） */
static inline double X_at_FortranIndex(const double *X, int idxF) {
    /* Fortran: X(1) 对应 C: X[0] */
    return X[idxF - 1];
}

/* Fortran 整数截断：idint 在这里可用 (int) 实现（对正数等价于 floor） */
static inline int idint_like(double a) {
    return (int)a;  /* trunc toward zero */
}

/* global_interp 的 C 版 */
void global_interp(const int ex[3],
                   const double *X, const double *Y, const double *Z,
                   const double *f,                 /* f(ex1,ex2,ex3) column-major */
                   double *f_int,
                   double x1, double y1, double z1,
                   int ORDN,
                   const double SoA[3],
                   int symmetry)
{
    enum { NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2 };

    int j, m;
    int imin, jmin, kmin;
    int cxB[3], cxT[3], cxI[3], cmin[3], cmax[3];
    double cx[3];
    double dX, dY, dZ, ddy;

    /* Fortran: imin=lbound(f,1) ... 通常是 1；这里按 1 处理 */
    imin = 1; jmin = 1; kmin = 1;

    dX = X_at_FortranIndex(X, imin + 1) - X_at_FortranIndex(X, imin);
    dY = X_at_FortranIndex(Y, jmin + 1) - X_at_FortranIndex(Y, jmin);
    dZ = X_at_FortranIndex(Z, kmin + 1) - X_at_FortranIndex(Z, kmin);

    /* x1a(j) = (j-1)*1.0  (j=1..ORDN) */
    double *x1a = (double*)malloc((size_t)ORDN * sizeof(double));
    double *ya  = (double*)malloc((size_t)ORDN * (size_t)ORDN * (size_t)ORDN * sizeof(double));
    if (!x1a || !ya) {
        fprintf(stderr, "global_interp: malloc failed\n");
        exit(1);
    }
    for (j = 0; j < ORDN; j++) x1a[j] = (double)j;

    /* cxI(m) = idint((p - P(1))/dP + 0.4) + 1  (Fortran 1-based) */
    cxI[0] = idint_like((x1 - X_at_FortranIndex(X, 1)) / dX + 0.4) + 1;
    cxI[1] = idint_like((y1 - X_at_FortranIndex(Y, 1)) / dY + 0.4) + 1;
    cxI[2] = idint_like((z1 - X_at_FortranIndex(Z, 1)) / dZ + 0.4) + 1;

    /* cxB = cxI - ORDN/2 + 1 ; cxT = cxB + ORDN - 1 */
    int half = ORDN / 2;  /* Fortran 整数除法 */
    for (m = 0; m < 3; m++) {
        cxB[m] = cxI[m] - half + 1;
        cxT[m] = cxB[m] + ORDN - 1;
    }

    /* cmin=1; cmax=ex */
    cmin[0] = cmin[1] = cmin[2] = 1;
    cmax[0] = ex[0];
    cmax[1] = ex[1];
    cmax[2] = ex[2];

    /* 对称边界时允许 cxB 为负/0（与 Fortran 一致） */
    if (symmetry == OCTANT && fabs(X_at_FortranIndex(X, 1)) < dX) cmin[0] = -half + 2;
    if (symmetry == OCTANT && fabs(X_at_FortranIndex(Y, 1)) < dY) cmin[1] = -half + 2;
    if (symmetry != NO_SYMM && fabs(X_at_FortranIndex(Z, 1)) < dZ) cmin[2] = -half + 2;

    /* 夹紧窗口 [cxB,cxT] 到 [cmin,cmax] */
    for (m = 0; m < 3; m++) {
        if (cxB[m] < cmin[m]) {
            cxB[m] = cmin[m];
            cxT[m] = cxB[m] + ORDN - 1;
        }
        if (cxT[m] > cmax[m]) {
            cxT[m] = cmax[m];
            cxB[m] = cxT[m] + 1 - ORDN;
        }
    }

    /*
      cx(m) 的计算：如果 cxB>0:
        cx = (p - P(cxB))/dP
      else:
        cx = (p + P(2 - cxB))/dP
      注意这里的 cxB 是 Fortran 1-based 语义下的整数，可能 <=0。
    */
    if (cxB[0] > 0) cx[0] = (x1 - X_at_FortranIndex(X, cxB[0])) / dX;
    else           cx[0] = (x1 + X_at_FortranIndex(X, 2 - cxB[0])) / dX;

    if (cxB[1] > 0) cx[1] = (y1 - X_at_FortranIndex(Y, cxB[1])) / dY;
    else           cx[1] = (y1 + X_at_FortranIndex(Y, 2 - cxB[1])) / dY;

    if (cxB[2] > 0) cx[2] = (z1 - X_at_FortranIndex(Z, cxB[2])) / dZ;
    else           cx[2] = (z1 + X_at_FortranIndex(Z, 2 - cxB[2])) / dZ;

    /* decide3d: 填充 ya(1:ORDN,1:ORDN,1:ORDN) */
    if (decide3d(ex, f, f, cxB, cxT, SoA, ya, ORDN, symmetry)) {
        printf("global_interp position: %g %g %g\n", x1, y1, z1);
        printf("data range: %g %g   %g %g   %g %g\n",
               X_at_FortranIndex(X, 1), X_at_FortranIndex(X, ex[0]),
               X_at_FortranIndex(Y, 1), X_at_FortranIndex(Y, ex[1]),
               X_at_FortranIndex(Z, 1), X_at_FortranIndex(Z, ex[2]));
        exit(1);
    }

    /* polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN) */
    polin3(x1a, x1a, x1a, ya, cx[0], cx[1], cx[2], f_int, &ddy, ORDN);

    free(x1a);
    free(ya);
}