#!/usr/bin/env bash
#SBATCH --job-name=ben1-pilot
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=ben1-pilot-%j.log
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 FO_JOBS=4 MPLBACKEND=Agg
unset DISPLAY
: "${LAB_REPO:?remote worktree path}"
: "${LAB_REVISION:?immutable pushed revision}"
: "${MEM_GUARD:?absolute path to cluster skill mem_guard.sh}"
cd "$LAB_REPO"
[[ $(git rev-parse HEAD) == "$LAB_REVISION" ]]
[[ -z $(git status --porcelain) ]]
runroot=${LAB_RUN_ROOT:-"$HOME/ben1-runs/$SLURM_JOB_ID"}
mkdir -p "$runroot"
cp "$MEM_GUARD" "$runroot/mem_guard.sh"
# Four enforced cores reserve 47.2 GiB. --mem is only documentation here.
export LAB_FLAGS='-O3 -fopenmp -finline-matmul-limit=0'
"$runroot/mem_guard.sh" 47.2 fo test --profile release --flag "$LAB_FLAGS" \
    test_lab_problem_definitions > "$runroot/build-test.log" 2>&1
# Fo's executable path is discovered from its supported exec interface once;
# the harness then times the binary without rebuilding or timing fo itself.
"$runroot/mem_guard.sh" 47.2 fo exec --profile release --flag "$LAB_FLAGS" \
    lab_scalar_run -- L1 1 /dev/null ichol 1 \
    > "$runroot/build-driver.log" 2>&1 && exit 1 || status=$?
# /dev/null intentionally cannot provide a mesh. The build must have completed.
# The binary is required below, and genuine execution failures are not rows.
LAB_EXE=$(find build -type f -name lab_scalar_run -perm -u+x | head -1)
[[ -n "$LAB_EXE" ]]
export LAB_EXE="$PWD/$LAB_EXE"
"$runroot/mem_guard.sh" 47.2 bash benchmark/lab/baselines.sh "$runroot/baselines.csv" \
    > "$runroot/runs.log" 2>&1
sha256sum "$runroot/baselines.csv" > "$runroot/SHA256SUMS"
