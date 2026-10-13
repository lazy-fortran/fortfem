#!/usr/bin/env bash
# Build the test programs with fpm and run them with bounded parallelism.
#
#   scripts/run_tests.sh [fast|slow|all] [-- extra fpm build flags]
#
#   fast  every test_*.f90 except *_slow (default; the CI selection)
#   slow  the test_*_slow.f90 programs
#   all   fast + slow (the run_* drivers are benchmarks and are not included)
#
# Environment: TEST_JOBS (concurrent programs), OMP_NUM_THREADS (threads per
# program, default 2), TEST_TIMEOUT (per program, default 120s), FPM_PROFILE
# (default release, as in CI).
#
# `fpm test` is not used to execute the programs: it starts all of them at once,
# so with several hundred OpenMP programs the machine is oversubscribed and
# programs hit their time limit although each runs in about a second alone.
set -u -o pipefail

group=${1:-fast}
[[ $# -gt 0 ]] && shift
[[ ${1:-} == -- ]] && shift

cd "$(dirname "$0")/.."

case $group in
fast) filter=(! -name '*_slow.f90') ;;
slow) filter=(-name '*_slow.f90') ;;
all) filter=() ;;
*) echo "usage: $0 [fast|slow|all]" >&2; exit 2 ;;
esac
mapfile -t names < <(find test -mindepth 1 -maxdepth 1 -type f -name 'test_*.f90' \
    "${filter[@]}" -printf '%f\n' | sed 's/\.f90$//' | sort)
if [[ ${#names[@]} -eq 0 ]]; then
    echo "no tests selected" >&2
    exit 2
fi

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-passive}
export TEST_TIMEOUT=${TEST_TIMEOUT:-120s}
# nproc itself honours OMP_NUM_THREADS, so ask it without.
cores=$(env -u OMP_NUM_THREADS nproc)
jobs=${TEST_JOBS:-$(( cores / OMP_NUM_THREADS ))}
(( jobs < 1 )) && jobs=1

# fpm builds the selected programs and, with `echo` as the runner, prints the
# path of each executable instead of running it.
paths=$(fpm test --profile "${FPM_PROFILE:-release}" --target "${names[@]}" "$@" \
    --runner echo | grep -E '/test/test_[A-Za-z0-9_]+$')
status=$?
if [[ $status -ne 0 || -z $paths ]]; then
    echo "fpm could not build the selected tests" >&2
    exit 1
fi

logdir=$(mktemp -d)
trap 'rm -rf "$logdir"' EXIT
export logdir

run_one() {
    local exe=$1 name
    name=$(basename "$exe")
    if timeout "$TEST_TIMEOUT" "$exe" >"$logdir/$name.log" 2>&1; then
        echo "PASS $name"
    else
        echo "FAIL $name (exit $?)"
        return 1
    fi
}
export -f run_one

echo "$paths" | xargs -P "$jobs" -I{} bash -c 'run_one {}' >"$logdir/summary"
failed=$(grep -c '^FAIL' "$logdir/summary")
passed=$(grep -c '^PASS' "$logdir/summary")
grep '^FAIL' "$logdir/summary" | sort | while read -r _ name _; do
    echo "--- $name"
    tail -n 15 "$logdir/$name.log"
done
echo "$passed passed, $failed failed of ${#names[@]} ($group, $jobs jobs x $OMP_NUM_THREADS threads)"
[[ $failed -eq 0 && $passed -eq ${#names[@]} ]]
