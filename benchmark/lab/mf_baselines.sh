#!/usr/bin/env bash
set -euo pipefail
: "${LAB_MF_EXE:?absolute path to lab_mf_run}"
: "${LAB_FLAGS:?exact compiler flags}"
out=${1:-mf.csv}
threads=${OMP_NUM_THREADS:-4}
for point in '2 1' '3 1' '4 1' '4 2'; do
    read -r p n <<< "$point"
    python3 benchmark/lab/L1.py --family gll-mf --p "$p" --n "$n" \
        --preconditioner pmg --threads "$threads" --executable "$LAB_MF_EXE" \
        --flags "$LAB_FLAGS" --output "$out"
done
failures=0
for ratio in 1 100 10000 1000000; do
    if ! python3 benchmark/lab/L3.py --family gll-mf --p 2 --n 1 --ratio "$ratio" \
        --preconditioner pmg --threads "$threads" --executable "$LAB_MF_EXE" \
        --flags "$LAB_FLAGS" --output "$out"; then
        # The failed row remains in the ledger; collection is not a pass gate.
        failures=$((failures + 1))
    fi
done
for n in 1 2; do
    python3 benchmark/lab/L5.py --family gll-mf --p 2 --n "$n" \
        --preconditioner pmg --threads "$threads" --executable "$LAB_MF_EXE" \
        --flags "$LAB_FLAGS" --output "$out"
    python3 benchmark/lab/L5.py --family gll-polar --p 2 --n "$n" \
        --preconditioner ichol_row --threads "$threads" --executable "$LAB_MF_EXE" \
        --flags "$LAB_FLAGS" --output "$out"
done
[[ $failures == 0 ]]
