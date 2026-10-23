#!/usr/bin/env bash
#SBATCH --job-name=seam1-l3
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --output=/home/%u/seam1-runs/%x-%j.log
# L3 aligned (mortar seam) vs unaligned rows on the exact unit-disk torus.
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 FO_JOBS=4 MPLBACKEND=Agg
unset DISPLAY
: "${LAB_REPO:?remote worktree path}"
: "${LAB_REVISION:?immutable pushed revision}"
: "${MEM_GUARD:?absolute path to cluster skill mem_guard.sh}"
cd "$LAB_REPO"
[[ $(git rev-parse HEAD) == "$LAB_REVISION" ]]
[[ -z $(git status --porcelain) ]]
runroot=${LAB_RUN_ROOT:-"$HOME/seam1-runs/$SLURM_JOB_ID"}
mkdir -p "$runroot"
cp "$MEM_GUARD" "$runroot/mem_guard.sh"
# Four enforced cores reserve 47.2 GiB. --mem is only documentation here.
export LAB_FLAGS='-O3 -funroll-loops -fopenmp -finline-matmul-limit=0'
# fo 0.3.2 links without -fopenmp, which fails on the cluster's offload-enabled
# GCC 12 (undefined GOMP_offload_unregister_ver); fpm passes the flags to the link.
"$runroot/mem_guard.sh" 47.2 fpm build --profile release --flag "$LAB_FLAGS" \
    > "$runroot/build-driver.log" 2>&1
mapfile -t executables < <(find build -path '*/app/lab_mf_run' -type f)
[[ ${#executables[@]} == 1 ]]
exe="$PWD/${executables[0]}"
points=${SEAM_L3_POINTS:-"2 4 6"}
# Ratio 1 is the isotropic control of the seam itself.
ratios=${SEAM_L3_RATIOS:-"1 1000 1000000"}
failures=0
for p in $points; do
    for ratio in $ratios; do
        for layout in disk aligned; do
            # A failed row stays in the ledger; collection is not a pass gate.
            "$runroot/mem_guard.sh" 47.2 python3 benchmark/lab/L3.py --family gll-mf \
                --p "$p" --n 1 --ratio "$ratio" --alignment "$layout" \
                --preconditioner pmg --threads 4 --executable "$exe" \
                --flags "$LAB_FLAGS" --scratch "$runroot" \
                --output "$runroot/seam-l3.csv" >> "$runroot/runs.log" 2>&1 \
                || failures=$((failures + 1))
        done
    done
done
sha256sum "$runroot/seam-l3.csv" > "$runroot/SHA256SUMS"
[[ $failures == 0 ]]
