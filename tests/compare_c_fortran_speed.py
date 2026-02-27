#!/usr/bin/env python3
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C_DIR = ROOT / "global_interp" / "c"
F_DIR = ROOT / "global_interp" / "fortran"
TEST_DIR = ROOT / "tests"
BUILD_DIR = ROOT / "build"


def run(cmd):
    return subprocess.run(cmd, check=True, text=True, capture_output=True)


def parse_elapsed(text: str) -> float:
    m = re.search(r"elapsed_sec=([0-9.]+)", text)
    if not m:
        raise RuntimeError(f"cannot parse elapsed from: {text}")
    return float(m.group(1))


def read_values(path: pathlib.Path):
    return [float(x.strip()) for x in path.read_text().splitlines() if x.strip()]


def main() -> int:
    BUILD_DIR.mkdir(exist_ok=True)
    gcc = shutil.which("gcc")
    gfortran = shutil.which("gfortran")
    if not gcc:
        print("ERROR: gcc not found", file=sys.stderr)
        return 2
    if not gfortran:
        print("SKIP: gfortran not found; cannot benchmark Fortran version in this environment")
        return 0

    c_exe = BUILD_DIR / "bench_c"
    f_exe = BUILD_DIR / "bench_fortran"

    run([
        gcc, "-O3", "-o", str(c_exe),
        str(TEST_DIR / "benchmark_c.c"),
        str(C_DIR / "global_interp.c"),
        str(C_DIR / "polint.c"),
        str(C_DIR / "polint3.c"),
        str(C_DIR / "decide3d.c"),
        "-lm",
    ])

    run([
        gfortran, "-O3", "-o", str(f_exe),
        str(TEST_DIR / "benchmark_fortran.f90"),
        str(F_DIR / "global_interp.f90"),
        str(F_DIR / "polint.f90"),
        str(F_DIR / "polint3.f90"),
        str(F_DIR / "decide3d.f90"),
    ])

    nquery = "15000"
    c_out = BUILD_DIR / "c_values.txt"
    f_out = BUILD_DIR / "fortran_values.txt"

    c_run = run([str(c_exe), nquery, str(c_out)])
    f_run = run([str(f_exe), nquery, str(f_out)])

    c_elapsed = parse_elapsed(c_run.stderr)
    f_elapsed = parse_elapsed(f_run.stdout)

    c_vals = read_values(c_out)
    f_vals = read_values(f_out)
    if len(c_vals) != len(f_vals):
        print(f"ERROR: value count mismatch C={len(c_vals)} Fortran={len(f_vals)}", file=sys.stderr)
        return 1

    max_abs = max(abs(a - b) for a, b in zip(c_vals, f_vals)) if c_vals else 0.0
    max_rel = max(abs(a - b) / max(1.0, abs(b)) for a, b in zip(c_vals, f_vals)) if c_vals else 0.0

    print(f"points={len(c_vals)}")
    print(f"max_abs_err={max_abs:.3e}")
    print(f"max_rel_err={max_rel:.3e}")
    print(f"c_elapsed_sec={c_elapsed:.6f}")
    print(f"fortran_elapsed_sec={f_elapsed:.6f}")
    print(f"speed_ratio_fortran_over_c={f_elapsed / c_elapsed:.3f}")

    if max_abs > 1e-10 and max_rel > 1e-10:
        print("FAILED: numerical mismatch beyond tolerance", file=sys.stderr)
        return 1
    print("PASS: C and Fortran outputs match within tolerance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
