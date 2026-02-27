#include <stdio.h>
#include <stdlib.h>

/*
  ex[0..2]  == Fortran ex(1:3)
  cxB/cxT   == Fortran cxB(1:3), cxT(1:3)  (可能 <=0)
  SoA[0..2] == Fortran SoA(1:3)
  f, fpi    == Fortran f(ex1,ex2,ex3) column-major (1-based in formulas)
  ya        == 连续内存，尺寸为 ORDN^3，对应 Fortran ya(cxB1:cxT1, cxB2:cxT2, cxB3:cxT3)
              但注意：我们用 offset 映射把 Fortran 的 i/j/k 坐标写进去。
*/

static inline int imax(int a, int b) { return a > b ? a : b; }
static inline int imin(int a, int b) { return a < b ? a : b; }

/* f(i,j,k): Fortran column-major, i/j/k are Fortran 1-based in [1..ex] */
#define F(i,j,k) f[((i)-1) + ex1 * (((j)-1) + ex2 * ((k)-1))]

/*
  ya(i,j,k): i in [cxB1..cxT1], j in [cxB2..cxT2], k in [cxB3..cxT3]
  我们把它映射到 C 的 0..ORDN-1 立方体：
    ii = i - cxB1
    jj = j - cxB2
    kk = k - cxB3
  并按 column-major 存储（与 Fortran 一致，方便直接喂给你的 polin3）
*/
#define YA(i,j,k) ya[((i)-cxB1) + ordn * (((j)-cxB2) + ordn * ((k)-cxB3))]

int decide3d(const int ex[3],
             const double *f,
             const double *fpi,   /* 这里未用，Fortran 也没用到 */
             const int cxB[3],
             const int cxT[3],
             const double SoA[3],
             double *ya,
             int ordn,
             int Symmetry)         /* Symmetry 在 decide3d 里也没直接用 */
{
    (void)fpi;
    (void)Symmetry;

    const int ex1 = ex[0], ex2 = ex[1], ex3 = ex[2];

    int fmin1[3], fmin2[3], fmax1[3], fmax2[3];
    int i, j, k, m;

    int gont = 0;

    /* 方便 YA 宏使用 */
    const int cxB1 = cxB[0], cxB2 = cxB[1], cxB3 = cxB[2];

    for (m = 0; m < 3; m++) {
        /* Fortran 的 “NaN 检查” 在整数上基本无意义，这里不额外处理 */

        fmin1[m] = imax(1, cxB[m]);
        fmax1[m] = cxT[m];

        fmin2[m] = cxB[m];
        fmax2[m] = imin(0, cxT[m]);

        /* if((fmin1<=fmax1) and (fmin1<1 or fmax1>ex)) gont=true */
        if ((fmin1[m] <= fmax1[m]) && (fmin1[m] < 1 || fmax1[m] > ex[m])) gont = 1;

        /* if((fmin2<=fmax2) and (2-fmax2<1 or 2-fmin2>ex)) gont=true */
        if ((fmin2[m] <= fmax2[m]) && (2 - fmax2[m] < 1 || 2 - fmin2[m] > ex[m])) gont = 1;
    }

    if (gont) {
        printf("error in decide3d\n");
        printf("cxB: %d %d %d   cxT: %d %d %d   ex: %d %d %d\n",
               cxB[0], cxB[1], cxB[2], cxT[0], cxT[1], cxT[2], ex[0], ex[1], ex[2]);
        printf("fmin1: %d %d %d  fmax1: %d %d %d\n",
               fmin1[0], fmin1[1], fmin1[2], fmax1[0], fmax1[1], fmax1[2]);
        printf("fmin2: %d %d %d  fmax2: %d %d %d\n",
               fmin2[0], fmin2[1], fmin2[2], fmax2[0], fmax2[1], fmax2[2]);
        return 1;
    }

    /* ---- 填充 ya：完全照 Fortran 两大块循环写 ---- */

    /* k in [fmin1(3)..fmax1(3)] */
    for (k = fmin1[2]; k <= fmax1[2]; k++) {

        /* j in [fmin1(2)..fmax1(2)] */
        for (j = fmin1[1]; j <= fmax1[1]; j++) {

            /* i in [fmin1(1)..fmax1(1)] : ya(i,j,k)=f(i,j,k) */
            for (i = fmin1[0]; i <= fmax1[0]; i++) {
                YA(i, j, k) = F(i, j, k);
            }

            /* i in [fmin2(1)..fmax2(1)] : ya(i,j,k)=f(2-i,j,k)*SoA(1) */
            for (i = fmin2[0]; i <= fmax2[0]; i++) {
                YA(i, j, k) = F(2 - i, j, k) * SoA[0];
            }
        }

        /* j in [fmin2(2)..fmax2(2)] */
        for (j = fmin2[1]; j <= fmax2[1]; j++) {

            /* i in [fmin1(1)..fmax1(1)] : ya(i,j,k)=f(i,2-j,k)*SoA(2) */
            for (i = fmin1[0]; i <= fmax1[0]; i++) {
                YA(i, j, k) = F(i, 2 - j, k) * SoA[1];
            }

            /* i in [fmin2(1)..fmax2(1)] : ya=f(2-i,2-j,k)*SoA(1)*SoA(2) */
            for (i = fmin2[0]; i <= fmax2[0]; i++) {
                YA(i, j, k) = F(2 - i, 2 - j, k) * SoA[0] * SoA[1];
            }
        }
    }

    /* k in [fmin2(3)..fmax2(3)] */
    for (k = fmin2[2]; k <= fmax2[2]; k++) {

        /* j in [fmin1(2)..fmax1(2)] */
        for (j = fmin1[1]; j <= fmax1[1]; j++) {

            /* i in [fmin1(1)..fmax1(1)] : ya=f(i,j,2-k)*SoA(3) */
            for (i = fmin1[0]; i <= fmax1[0]; i++) {
                YA(i, j, k) = F(i, j, 2 - k) * SoA[2];
            }

            /* i in [fmin2(1)..fmax2(1)] : ya=f(2-i,j,2-k)*SoA(1)*SoA(3) */
            for (i = fmin2[0]; i <= fmax2[0]; i++) {
                YA(i, j, k) = F(2 - i, j, 2 - k) * SoA[0] * SoA[2];
            }
        }

        /* j in [fmin2(2)..fmax2(2)] */
        for (j = fmin2[1]; j <= fmax2[1]; j++) {

            /* i in [fmin1(1)..fmax1(1)] : ya=f(i,2-j,2-k)*SoA(2)*SoA(3) */
            for (i = fmin1[0]; i <= fmax1[0]; i++) {
                YA(i, j, k) = F(i, 2 - j, 2 - k) * SoA[1] * SoA[2];
            }

            /* i in [fmin2(1)..fmax2(1)] : ya=f(2-i,2-j,2-k)*SoA1*SoA2*SoA3 */
            for (i = fmin2[0]; i <= fmax2[0]; i++) {
                YA(i, j, k) = F(2 - i, 2 - j, 2 - k) * SoA[0] * SoA[1] * SoA[2];
            }
        }
    }

    return 0;
}

#undef F
#undef YA