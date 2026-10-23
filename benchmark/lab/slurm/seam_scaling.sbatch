#!/usr/bin/env bash
#SBATCH --job-name=seam1-scal
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --output=/home/%u/seam1-runs/%x-%j.log
# Strong scaling of the aligned L3 mortar seam setup, isotropic ratio 1, n=1.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg OMP_PROC_BIND=close OMP_PLACES=cores
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
# Each enforced core reserves 11.8 GiB. --mem is only documentation here.
guard=$(awk -v c="${SLURM_CPUS_PER_TASK:-16}" 'BEGIN {printf "%.1f", 11.8*c}')
# The cluster's GCC 12 is offload-enabled, see seam_l3.sbatch.
export LAB_FLAGS='-O3 -funroll-loops -fopenmp -finline-matmul-limit=0 -foffload=disable'
fpm_args=(--profile release --flag "$LAB_FLAGS" --link-flag "-fopenmp -foffload=disable")
FO_JOBS=16 "$runroot/mem_guard.sh" "$guard" fpm run "${fpm_args[@]}" \
    --target lab_mf_run --runner echo > "$runroot/build-driver.log" 2>&1
exe=$(grep -E '(^|/)lab_mf_run$' "$runroot/build-driver.log" | tail -n 1)
[[ $exe == /* ]] || exe="$PWD/$exe"
[[ -x $exe ]]
sha256sum "$exe" > "$runroot/executable.sha256"
out="$runroot/seam-scaling.txt"
echo "# p threads seam_setup_s wall_s peak_rss_kb rel_l2 iterations host" > "$out"
for p in ${SEAM_SCALING_POINTS:-"4 6"}; do
    for t in ${SEAM_SCALING_THREADS:-"1 2 4 8 16"}; do
        log="$runroot/p${p}-t${t}.log"
        # The nodes have no /usr/bin/time: wall time and the child high-water
        # mark come from wait4, as in run.py.
        OMP_NUM_THREADS=$t python3 -c '
import os, subprocess, sys, time
t0 = time.perf_counter()
proc = subprocess.Popen(sys.argv[1:])
_, status, usage = os.wait4(proc.pid, 0)
print("TIME", time.perf_counter() - t0, usage.ru_maxrss, flush=True)
sys.exit(os.waitstatus_to_exitcode(status))
' "$runroot/mem_guard.sh" "$guard" "$exe" L3 "$p" 1 pmg 1 aligned > "$log" 2>&1
        seam=$(awk '/^SEAM/ {print $3}' "$log")
        read -r wall rss < <(awk '/^TIME/ {print $2, $3}' "$log")
        read -r err its < <(awk '/^LAB/ {print $3, $9}' "$log")
        echo "$p $t $seam $wall $rss $err $its $(hostname)" >> "$out"
    done
done
sha256sum "$out" > "$runroot/SHA256SUMS"
