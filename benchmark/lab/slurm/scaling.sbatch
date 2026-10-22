#!/usr/bin/env bash
#SBATCH --job-name=ben1-scale
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --array=0-6%1
#SBATCH --output=ben1-scale-%A-%a.log
set -euo pipefail
# Submit only after the pilot's observed peak is known. Do not shrink a
# memory-driven core reservation below ceil(1.5*peak_GiB/11.8).
# This template reserves 64 cores for each measurement, including smaller
# thread counts, to hold placement fixed. It is an explicit one-node sweep.
threads=$((1 << SLURM_ARRAY_TASK_ID))
export OMP_NUM_THREADS=$threads OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
unset DISPLAY
: "${LAB_REPO:?}" "${LAB_EXE:?}" "${LAB_FLAGS:?}" "${MEM_GUARD:?}"
: "${LAB_FAMILY:?tet, bspline or gll-mf}"
: "${LAB_P:?}" "${LAB_N:?pilot-qualified base resolution}"
: "${LAB_MODE:?strong or weak}"
cd "$LAB_REPO"
n=$LAB_N
case "$LAB_MODE" in
strong) ;;
weak)
    # Fixed p, 3D weak scaling: actual DOFs/thread reported, because integer
    # mesh rounding and constraints prevent exactly constant DOFs/thread.
    n=$(python3 -c 'import math,sys; print(max(1,round(int(sys.argv[1])*int(sys.argv[2])**(1/3))))' "$n" "$threads") ;;
*) exit 2 ;;
esac
"$MEM_GUARD" 755.2 python3 benchmark/lab/L1.py --family "$LAB_FAMILY" \
    --p "$LAB_P" --n "$n" --preconditioner ichol --threads "$threads" \
    --executable "$LAB_EXE" --flags "$LAB_FLAGS" \
    --output "${LAB_MODE}-${LAB_FAMILY}-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.csv"
