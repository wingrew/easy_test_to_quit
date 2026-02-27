# Numerical validation & performance comparison

## 1) Validate C numerical correctness

This check compares the C implementation against a Python translation of the Fortran logic
(`global_interp` + `decide3d` + `polin3/polint`) on deterministic test points.

```bash
python3 tests/validate_c_numerics.py
```

## 2) Compare C vs Fortran speed and outputs

This benchmark compiles and runs both implementations on the same dataset and random query sequence,
then reports max error and runtime ratio.

```bash
python3 tests/compare_c_fortran_speed.py
```

If `gfortran` is unavailable, the script reports a skip message.
