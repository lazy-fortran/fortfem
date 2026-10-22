#!/usr/bin/env bash
#SBATCH --job-name=ben1-scale
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --output=ben1-scale-%j.log
set -euo pipefail
# Submit separately with --cpus-per-task=max(threads,ceil(1.5*pilot_peak_GiB/11.8)).
# Set LAB_THREADS in {1,2,4,8,16,32,64}; do not run with fewer allocated cores.
: "${LAB_THREADS:?}" "${LAB_PILOT_PEAK_GIB:?measured completed pilot peak}"
threads=$LAB_THREADS
[[ $threads -le $SLURM_CPUS_PER_TASK ]]
required=$(python3 -c 'import math,sys; print(math.ceil(1.5*float(sys.argv[1])/11.8))' "$LAB_PILOT_PEAK_GIB")
[[ $required -le $SLURM_CPUS_PER_TASK ]]
limit=$(python3 -c 'import sys; print(int(sys.argv[1])*11.8)' "$SLURM_CPUS_PER_TASK")
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
"$MEM_GUARD" "$limit" python3 benchmark/lab/L1.py --family "$LAB_FAMILY" \
    --p "$LAB_P" --n "$n" --preconditioner ichol --threads "$threads" \
    --executable "$LAB_EXE" --flags "$LAB_FLAGS" \
    --output "${LAB_MODE}-${LAB_FAMILY}-${SLURM_JOB_ID}-${threads}.csv"
