#!/usr/bin/env python3
import ctypes
import math
import pathlib
import subprocess
import sys
from typing import Tuple

ROOT = pathlib.Path(__file__).resolve().parents[1]
C_DIR = ROOT / "global_interp" / "c"
BUILD_DIR = ROOT / "build"


def build_c_shared() -> pathlib.Path:
    BUILD_DIR.mkdir(exist_ok=True)
    so_path = BUILD_DIR / "libglobal_interp_c.so"
    cmd = [
        "gcc", "-O3", "-shared", "-fPIC", "-o", str(so_path),
        str(C_DIR / "global_interp.c"),
        str(C_DIR / "polint.c"),
        str(C_DIR / "polint3.c"),
        str(C_DIR / "decide3d.c"),
        "-lm",
    ]
    subprocess.run(cmd, check=True)
    return so_path


def col_major_idx(i: int, j: int, k: int, ex1: int, ex2: int) -> int:
    return i + ex1 * (j + ex2 * k)


def decide3d_py(ex, f, cxB, cxT, soa, ordn):
    ex1, ex2, ex3 = ex
    cxB1, cxB2, cxB3 = cxB

    def ya_idx(i, j, k):
        return (i - cxB1) + ordn * ((j - cxB2) + ordn * (k - cxB3))

    def F(i, j, k):
        return f[col_major_idx(i - 1, j - 1, k - 1, ex1, ex2)]

    fmin1 = [max(1, cxB[m]) for m in range(3)]
    fmax1 = [cxT[m] for m in range(3)]
    fmin2 = [cxB[m] for m in range(3)]
    fmax2 = [min(0, cxT[m]) for m in range(3)]

    for m in range(3):
        if (fmin1[m] <= fmax1[m]) and (fmin1[m] < 1 or fmax1[m] > ex[m]):
            raise RuntimeError("decide3d_py out of range")
        if (fmin2[m] <= fmax2[m]) and (2 - fmax2[m] < 1 or 2 - fmin2[m] > ex[m]):
            raise RuntimeError("decide3d_py mirror out of range")

    ya = [0.0] * (ordn * ordn * ordn)

    for k in range(fmin1[2], fmax1[2] + 1):
        for j in range(fmin1[1], fmax1[1] + 1):
            for i in range(fmin1[0], fmax1[0] + 1):
                ya[ya_idx(i, j, k)] = F(i, j, k)
            for i in range(fmin2[0], fmax2[0] + 1):
                ya[ya_idx(i, j, k)] = F(2 - i, j, k) * soa[0]
        for j in range(fmin2[1], fmax2[1] + 1):
            for i in range(fmin1[0], fmax1[0] + 1):
                ya[ya_idx(i, j, k)] = F(i, 2 - j, k) * soa[1]
            for i in range(fmin2[0], fmax2[0] + 1):
                ya[ya_idx(i, j, k)] = F(2 - i, 2 - j, k) * soa[0] * soa[1]

    for k in range(fmin2[2], fmax2[2] + 1):
        for j in range(fmin1[1], fmax1[1] + 1):
            for i in range(fmin1[0], fmax1[0] + 1):
                ya[ya_idx(i, j, k)] = F(i, j, 2 - k) * soa[2]
            for i in range(fmin2[0], fmax2[0] + 1):
                ya[ya_idx(i, j, k)] = F(2 - i, j, 2 - k) * soa[0] * soa[2]
        for j in range(fmin2[1], fmax2[1] + 1):
            for i in range(fmin1[0], fmax1[0] + 1):
                ya[ya_idx(i, j, k)] = F(i, 2 - j, 2 - k) * soa[1] * soa[2]
            for i in range(fmin2[0], fmax2[0] + 1):
                ya[ya_idx(i, j, k)] = F(2 - i, 2 - j, 2 - k) * soa[0] * soa[1] * soa[2]

    return ya


def polint_py(xa, ya, x):
    ordn = len(xa)
    c = ya.copy()
    d = ya.copy()
    ho = [xa_i - x for xa_i in xa]
    ns_f = 1
    dif = abs(x - xa[0])
    for i in range(1, ordn):
        dift = abs(x - xa[i])
        if dift < dif:
            ns_f = i + 1
            dif = dift
    y = ya[ns_f - 1]
    ns_f -= 1
    for m in range(1, ordn):
        n_m = ordn - m
        for i in range(n_m):
            hp = ho[i]
            h = ho[i + m]
            den = hp - h
            den = (c[i + 1] - d[i]) / den
            d[i] = h * den
            c[i] = hp * den
        if 2 * ns_f < n_m:
            dy = c[ns_f]
        else:
            dy = d[ns_f - 1]
            ns_f -= 1
        y += dy
    return y


def polin3_py(x1a, x2a, x3a, ya, x1, x2, x3, ordn):
    def ya3(i, j, k):
        return ya[i + ordn * (j + ordn * k)]

    yatmp = [0.0] * (ordn * ordn)
    ymtmp = [0.0] * ordn
    for k in range(ordn):
        for j in range(ordn):
            vec = [ya3(i, j, k) for i in range(ordn)]
            yatmp[j + ordn * k] = polint_py(x1a, vec, x1)
    for k in range(ordn):
        vec = [yatmp[j + ordn * k] for j in range(ordn)]
        ymtmp[k] = polint_py(x2a, vec, x2)
    return polint_py(x3a, ymtmp, x3)


def global_interp_py(ex, X, Y, Z, f, x1, y1, z1, ordn, soa, symmetry):
    no_symm, octant = 0, 2
    dx = X[1] - X[0]
    dy = Y[1] - Y[0]
    dz = Z[1] - Z[0]

    def idint_like(a):
        return int(a)

    cxI = [
        idint_like((x1 - X[0]) / dx + 0.4) + 1,
        idint_like((y1 - Y[0]) / dy + 0.4) + 1,
        idint_like((z1 - Z[0]) / dz + 0.4) + 1,
    ]
    half = ordn // 2
    cxB = [cxI[m] - half + 1 for m in range(3)]
    cxT = [cxB[m] + ordn - 1 for m in range(3)]

    cmin = [1, 1, 1]
    cmax = [ex[0], ex[1], ex[2]]
    if symmetry == octant and abs(X[0]) < dx:
        cmin[0] = -half + 2
    if symmetry == octant and abs(Y[0]) < dy:
        cmin[1] = -half + 2
    if symmetry != no_symm and abs(Z[0]) < dz:
        cmin[2] = -half + 2

    for m in range(3):
        if cxB[m] < cmin[m]:
            cxB[m] = cmin[m]
            cxT[m] = cxB[m] + ordn - 1
        if cxT[m] > cmax[m]:
            cxT[m] = cmax[m]
            cxB[m] = cxT[m] + 1 - ordn

    cx = [0.0, 0.0, 0.0]
    cx[0] = (x1 - X[cxB[0] - 1]) / dx if cxB[0] > 0 else (x1 + X[(2 - cxB[0]) - 1]) / dx
    cx[1] = (y1 - Y[cxB[1] - 1]) / dy if cxB[1] > 0 else (y1 + Y[(2 - cxB[1]) - 1]) / dy
    cx[2] = (z1 - Z[cxB[2] - 1]) / dz if cxB[2] > 0 else (z1 + Z[(2 - cxB[2]) - 1]) / dz

    x1a = [float(i) for i in range(ordn)]
    ya = decide3d_py(ex, f, cxB, cxT, soa, ordn)
    return polin3_py(x1a, x1a, x1a, ya, cx[0], cx[1], cx[2], ordn)


def make_data() -> Tuple[Tuple[int, int, int], list, list, list, list]:
    ex = (16, 14, 12)
    dx, dy, dz = 0.1, 0.2, 0.15
    X = [i * dx for i in range(ex[0])]
    Y = [j * dy for j in range(ex[1])]
    Z = [k * dz for k in range(ex[2])]
    f = [0.0] * (ex[0] * ex[1] * ex[2])
    for k in range(ex[2]):
        for j in range(ex[1]):
            for i in range(ex[0]):
                x, y, z = X[i], Y[j], Z[k]
                f[col_major_idx(i, j, k, ex[0], ex[1])] = math.sin(x) + 0.25 * math.cos(1.7 * y) + z * z + 0.05 * x * y
    return ex, X, Y, Z, f


def main() -> int:
    so = build_c_shared()
    lib = ctypes.CDLL(str(so))
    fn = lib.global_interp
    fn.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]

    ex, X, Y, Z, f = make_data()
    ordn = 4
    soa = [1.0, -1.0, 1.0]
    symmetry = 2

    ex_a = (ctypes.c_int * 3)(*ex)
    X_a = (ctypes.c_double * len(X))(*X)
    Y_a = (ctypes.c_double * len(Y))(*Y)
    Z_a = (ctypes.c_double * len(Z))(*Z)
    f_a = (ctypes.c_double * len(f))(*f)
    soa_a = (ctypes.c_double * 3)(*soa)

    max_abs = 0.0
    max_rel = 0.0
    bad = 0

    for q in range(200):
        xq = (q * 37 % 997) / 997.0 * X[-1]
        yq = (q * 53 % 991) / 991.0 * Y[-1]
        zq = (q * 71 % 983) / 983.0 * Z[-1]

        out_c = ctypes.c_double(0.0)
        fn(ex_a, X_a, Y_a, Z_a, f_a, ctypes.byref(out_c), xq, yq, zq, ordn, soa_a, symmetry)
        out_py = global_interp_py(ex, X, Y, Z, f, xq, yq, zq, ordn, soa, symmetry)

        abs_err = abs(out_c.value - out_py)
        rel_err = abs_err / max(1.0, abs(out_py))
        max_abs = max(max_abs, abs_err)
        max_rel = max(max_rel, rel_err)
        if abs_err > 1e-10 and rel_err > 1e-10:
            bad += 1

    print(f"max_abs_err={max_abs:.3e}")
    print(f"max_rel_err={max_rel:.3e}")
    if bad:
        print(f"FAILED: {bad} points exceed tolerance", file=sys.stderr)
        return 1
    print("PASS: C implementation matches Python-translated Fortran logic")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
