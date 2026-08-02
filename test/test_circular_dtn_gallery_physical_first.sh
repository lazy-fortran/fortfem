#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/circular_dtn_modes"
log_file="${TMPDIR:-/tmp}/fortfem-circular-dtn-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" \
        fo exec --no-build circular_dtn_modes >"$log_file"
)

test -s "$output_directory/circular_dtn_field_2d.png"
test -s "$output_directory/circular_dtn_trace_1d.png"
test -s "$output_directory/circular_dtn_response_1d.png"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/circular_dtn_field_2d.png" | \
    cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "circular_dtn_modes" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = circular_dtn_field_2d.png

# The ninth boundary harmonic is deliberately truncated by the retained mode
# set.  Check that independent analytical value in addition to plot staging.
python3 - "$log_file" <<'PY'
import math
import pathlib
import re
import sys

log = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
match = re.search(r"discarded relative trace norm:\s*([0-9.eE+-]+)", log)
assert match, "DtN gallery did not report its discarded-mode diagnostic"
value = float(match.group(1))
expected = 0.1 / math.sqrt(1.01)
assert abs(value - expected) < 2.0e-5
PY

echo "circular DtN gallery renders its physical disk field before diagnostics"
