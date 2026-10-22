#!/usr/bin/env bash
#SBATCH --job-name=ben1-pilot
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=/home/%u/ben1-runs/%x-%j.log
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
export LAB_FLAGS='-O3 -funroll-loops -fopenmp -finline-matmul-limit=0'
"$runroot/mem_guard.sh" 47.2 fo test --profile release --flag "$LAB_FLAGS" \
    test_lab_problem_definitions > "$runroot/build-test.log" 2>&1
"$runroot/mem_guard.sh" 47.2 fo build --profile release --flag "$LAB_FLAGS" \
    > "$runroot/build-driver.log" 2>&1
# A fresh, dedicated pilot worktree has one profile. Ambiguity is an error,
# never a reason to choose an old or newest-mtime executable.
mapfile -t executables < <(find build/fo/profiles -path '*/bin/lab_scalar_run' -type f)
[[ ${#executables[@]} == 1 ]]
export LAB_EXE="$PWD/${executables[0]}"
"$runroot/mem_guard.sh" 47.2 bash benchmark/lab/baselines.sh "$runroot/baselines.csv" \
    > "$runroot/runs.log" 2>&1
sha256sum "$runroot/baselines.csv" > "$runroot/SHA256SUMS"
