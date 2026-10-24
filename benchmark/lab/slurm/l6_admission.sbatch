#!/usr/bin/env bash
#SBATCH --job-name=ben2-l6-admit
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=/home/%u/ben1-runs/%x-%j.log
set -euo pipefail
# Reads the retained N1 records read-only and writes only under the run root.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
unset DISPLAY
: "${LAB_REPO:?remote worktree path}" "${LAB_REVISION:?immutable pushed revision}"
: "${MEM_GUARD:?absolute path to cluster skill mem_guard.sh}"
: "${N1_ROOT:?retained N1 tree, read-only}" "${N1_PYTHON:?N1 SIMSOPT python}"
cd "$LAB_REPO"
[[ $(git rev-parse HEAD) == "$LAB_REVISION" ]]
[[ -z $(git status --porcelain) ]]
runroot=${LAB_RUN_ROOT:-"$HOME/ben1-runs/$SLURM_JOB_ID"}
mkdir -p "$runroot"
cp "$MEM_GUARD" "$runroot/mem_guard.sh"
# Four enforced cores reserve 47.2 GiB.
"$runroot/mem_guard.sh" 47.2 "$N1_PYTHON" -I benchmark/lab/import_n1.py "$N1_ROOT" \
    --admit --output "$runroot/l6-n1-admission.csv" > "$runroot/admission.log" 2>&1
sha256sum "$runroot/l6-n1-admission.csv" > "$runroot/SHA256SUMS"
