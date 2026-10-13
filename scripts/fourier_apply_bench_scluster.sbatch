#!/bin/bash
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 01:30:00
#SBATCH --mem=24G
#SBATCH -o /home/ert/fortfem-fft-%j.log
# Matrix-free Fourier apply: test, then memory/time vs assembled form, N = 4..16.
# Usage: cd $HOME/fortfem-fft-src; sbatch scripts/fourier_apply_bench_scluster.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR}"
E=$HOME/fs1/env
export PATH="$HOME/fsym-bin:$E/bin:$HOME/.local/bin:$PATH" LD_LIBRARY_PATH="$E/lib:${LD_LIBRARY_PATH:-}" LIBRARY_PATH="$E/lib:$HOME/fsym-env/lib" CPATH="$HOME/fsym-gmpinc:$E/include"
export FPM_FC=gfortran FPM_CC=gcc FPM_CXX=g++ OMP_PLACES=cores OMP_PROC_BIND=spread OMP_DYNAMIC=false
OMP_NUM_THREADS=4 fpm test test_meridian_fourier_apply --profile release > test.log 2>&1
echo "test exit $?" >> test.log
fpm build --profile release > build.log 2>&1
: > bench.log
for t in 1 4 16; do
  OMP_NUM_THREADS=$t fpm run --profile release fortfem_fourier_apply_bench -- 16 4 8 12 16 >> bench.log 2>&1
done
