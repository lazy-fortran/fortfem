#!/usr/bin/env bash
# Run only on a permitted remote compute host. Build is excluded from timings.
set -euo pipefail
: "${LAB_EXE:?absolute path to lab_scalar_run}"
: "${LAB_FLAGS:?exact flags used to build LAB_EXE}"
out=${1:-lab.csv}
threads=${OMP_NUM_THREADS:-4}
for p in 1 2 3; do
    for n in 1 2; do
        python3 benchmark/lab/L1.py --family tet --p "$p" --n "$n" \
            --preconditioner ichol_row --threads "$threads" \
            --executable "$LAB_EXE" --flags "$LAB_FLAGS" --output "$out"
    done
done
for n in 8 16 32 64; do
    for ratio in 1 100 10000 1000000; do
        python3 benchmark/lab/L3a.py --family tet --p 1 --n "$n" \
            --ratio "$ratio" --preconditioner ldlt --threads "$threads" \
            --executable "$LAB_EXE" --flags "$LAB_FLAGS" --output "$out"
    done
done
for ratio in 1 100 10000 1000000; do
    python3 benchmark/lab/L3.py --family tet --p 2 --n 1 \
        --ratio "$ratio" --preconditioner ldlt --threads "$threads" \
        --executable "$LAB_EXE" --flags "$LAB_FLAGS" --output "$out"
done

# One larger assembled point, still below a few 1e5 unknowns.
python3 benchmark/lab/L1.py --family tet --p 3 --n 3 \
    --preconditioner ichol_row --threads "$threads" \
    --executable "$LAB_EXE" --flags "$LAB_FLAGS" --output "$out"
if [[ -n ${LAB_SPLINE_EXE:-} ]]; then
    for p in 2 3 4; do
        for n in 1 2; do
            python3 benchmark/lab/L1.py --family bspline --p "$p" --n "$n" \
                --preconditioner ichol_row --threads "$threads" \
                --executable "$LAB_SPLINE_EXE" --flags "$LAB_FLAGS" --output "$out"
        done
    done
fi
