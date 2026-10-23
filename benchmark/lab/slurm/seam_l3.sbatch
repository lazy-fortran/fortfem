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
# The cluster's GCC 12 is offload-enabled: without -foffload=disable its link
# adds an offload stub that leaves GOMP_offload_unregister_ver undefined (fo
# 0.3.2 and fpm alike). fpm resolves the binary for exactly these flags.
export LAB_FLAGS='-O3 -funroll-loops -fopenmp -finline-matmul-limit=0 -foffload=disable'
fpm_args=(--profile release --flag "$LAB_FLAGS" --link-flag "-fopenmp -foffload=disable")
"$runroot/mem_guard.sh" 47.2 fpm run "${fpm_args[@]}" --target lab_mf_run \
    --runner echo > "$runroot/build-driver.log" 2>&1
# mem_guard appends its peak line; take the path the runner echoed.
exe=$(grep -E '(^|/)lab_mf_run$' "$runroot/build-driver.log" | tail -n 1)
[[ $exe == /* ]] || exe="$PWD/$exe"
[[ -x $exe ]]
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
