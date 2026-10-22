#!/bin/bash
# NEQ1 lab Grad-Shafranov timing on one exclusive sCluster node.
# Submit from anywhere:  sbatch --export=ALL,MODE=pilot|scan lab.sbatch
# Logs go to $RUN/logs (outside the worktree); REPO must be a clean checkout.
#SBATCH --job-name=neq1-lab
#SBATCH --partition=compute
#SBATCH --exclusive
#SBATCH --cpus-per-task=96
#SBATCH --time=02:00:00
#SBATCH --output=/home/ert/neq1-20261022/logs/lab-%j.log
set -Eeuo pipefail
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
unset DISPLAY
RUN=${RUN:-/home/ert/neq1-20261022}
REPO=${REPO:-/home/ert/wt/fortfem/neq1-20261022}
PY=${PY:-/home/ert/k6bench-20261021/venv/bin/python}
MODE=${MODE:-pilot}
BUILD=$RUN/lab/build
CASE=$RUN/vmec/case
OUT=$RUN/lab/$MODE-$SLURM_JOB_ID
mkdir -p "$OUT"
echo "HOST=$(hostname) JOB=$SLURM_JOB_ID MODE=$MODE REVISION=$(git -C "$REPO" rev-parse HEAD)"
[[ -z $(git -C "$REPO" status --porcelain) ]] || { echo "dirty worktree"; exit 2; }
grep -m1 'model name' /proc/cpuinfo; gfortran --version | head -1
sha256sum /lib/x86_64-linux-gnu/libc.so.6
FLAGS='-O3 -march=native -funroll-loops'
"$RUN/vmec/mem_guard.sh" 64 cmake -S "$REPO" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_FLAGS="$FLAGS" -DFORTFEM_BUILD_TESTING=ON \
    -DFORTFEM_BUILD_BENCHMARKS=ON -DFORTFEM_NATIVE_SPARSE_ONLY=ON > "$OUT/configure.log"
"$RUN/vmec/mem_guard.sh" 64 cmake --build "$BUILD" -j16 --target fortfem_lab_neq1 \
    test_lab_newton_krylov test_lab_grad_shafranov test_lab_tensor_mesh \
    test_lab_tensor_apply test_lab_tensor_solver test_lab_gll_line > "$OUT/build.log"
sha256sum "$BUILD/fortfem_lab_neq1" | tee "$OUT/binary.sha256"
echo "FLAGS=$FLAGS"
OMP_NUM_THREADS=4 ctest --test-dir "$BUILD" --output-on-failure \
    -R '^test_lab_(gll_line|tensor_mesh|tensor_apply|tensor_solver|newton_krylov|grad_shafranov)$'

# point CASEDIR n p threads tag [rtol]
point() {
    local casedir=$1 n=$2 p=$3 t=$4 tag=$5 rtol=${6:-1e-10}
    local dir=$OUT/$tag/n$n-p$p-t$t cpus
    mkdir -p "$dir"
    cpus=$(seq -s, 0 $((t - 1)))
    OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores "$RUN/vmec/mem_guard.sh" 64 \
        /usr/bin/time -v taskset -c "$cpus" "$BUILD/fortfem_lab_neq1" \
        "$casedir/case.nml" "$n" "$p" "$casedir/samples.txt" "$dir/" "$rtol" 128 \
        > "$dir/solve.log" 2> "$dir/time.txt" || echo "FAILED $dir"
    if [[ -f $casedir/reference_samples.txt && -f $dir/samples_out.txt ]]; then
        "$PY" -I "$REPO/benchmark/lab/neq1/l6_metric.py" "$casedir/reference_samples.txt" \
            "$dir/samples_out.txt" > "$dir/l6.json" || true
    fi
    echo "$tag n=$n p=$p t=$t $(grep -o 't_solve=.*' "$dir/solve.log") $(cat "$dir/l6.json" 2>/dev/null | tr -d '\n' | cut -c1-160)"
}

if [[ $MODE == pilot ]]; then
    point "$CASE" 2 6 1 solovev
    point "$CASE" 4 6 8 solovev
    touch "$OUT/complete"
    exit 0
fi

# Case A (exact Solov'ev): repeat each point three times (rep tag) for timing noise.
for rep in 1 2 3; do
    for t in 1 8; do
        for n in 1 2 3 4 6 8; do
            for p in 2 3 4 5 6 7 8 10; do
                (( n * p > 64 )) && continue
                point "$CASE" "$n" "$p" "$t" "solovev-r$rep"
            done
        done
    done
done
touch "$OUT/complete"
