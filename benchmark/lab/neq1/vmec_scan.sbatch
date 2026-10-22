#!/bin/bash
#SBATCH --job-name=neq1-vmec
#SBATCH --partition=compute
#SBATCH --exclusive
#SBATCH --cpus-per-task=96
#SBATCH --nodelist=node12
#SBATCH --time=06:00:00
#SBATCH --output=/home/ert/neq1-20261022/vmec/logs/%x-%j.log
# VMEC++ matched-accuracy scan on the exact Solov'ev case.
# MODE=pilot (one point) or MODE=scan (default):  sbatch --export=ALL,MODE=pilot vmec_scan.sbatch
# Points run serially on an exclusive node, each pinned to the first `threads` cores,
# under mem_guard (64 GiB) and /usr/bin/time -v.  Resumable: finished points are skipped.
set -Eeuo pipefail
ROOT=/home/ert/neq1-20261022/vmec
cd "$ROOT"
PY=/home/ert/k6bench-20261021/venv/bin/python
MODE=${MODE:-scan}
export MEM_GUARD_INTERVAL=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
unset DISPLAY
OUTROOT=$ROOT/$MODE
mkdir -p "$OUTROOT"
{ hostname; date -Is; uname -a; lscpu | head -20; echo "SLURM_JOB_ID=$SLURM_JOB_ID";
  $PY -c 'import vmecpp,simsopt,numpy,scipy;print("vmecpp",vmecpp.__version__ if hasattr(vmecpp,"__version__") else "?");print("numpy",numpy.__version__,"scipy",scipy.__version__)';
  sha256sum py/*.py case/profiles.txt case/boundary.txt case/samples.txt case/reference_samples.txt
} > "$OUTROOT/machine.txt" 2>&1

runpoint() {
 local mpol=$1 ns=$2 ftol=$3 threads=$4
 local out="$OUTROOT/m${mpol}_ns${ns}_f${ftol}_t${threads}"
 if [ -s "$out/l6.json" ]; then echo "skip $out"; return; fi
 mkdir -p "$out"
 local cpus
 cpus=$(/usr/bin/python3 -c 'import os,sys;print(",".join(map(str,sorted(os.sched_getaffinity(0))[:int(sys.argv[1])])))' "$threads")
 export OMP_NUM_THREADS=$threads
 set +e
 ./mem_guard.sh 64 /usr/bin/time -v -o "$out/time.txt" timeout 1200 taskset -c "$cpus" \
   $PY py/vmec_case.py --case case --mpol "$mpol" --ns "$ns" --ftol "$ftol" \
   --threads "$threads" --out "$out" > "$out/solve.log" 2> "$out/stderr.log"
 local rc=$?
 echo "$rc" > "$out/exit.txt"
 set -e
 echo "$out exit=$rc"
}

if [ "$MODE" = pilot ]; then
 runpoint 8 51 1e-12 1
else
 # high-resolution convention-check point first (psi_axis, e_flux), then the full grid
 runpoint 16 201 1e-14 8
 for ftol in 1e-11 1e-14; do
  for ns in 25 51 101 201; do
   for mpol in 4 6 8 10 12 16; do
    for t in 1 8; do runpoint "$mpol" "$ns" "$ftol" "$t"; done
   done
  done
 done
fi
touch "$OUTROOT/complete"
