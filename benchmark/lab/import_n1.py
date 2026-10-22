#!/usr/bin/env python3
"""Read the existing N1 sweep ledger; never rerun or copy equilibrium cases."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import sys

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('root', type=Path)
p.add_argument('--output', required=True)
a = p.parse_args()
root = a.root.resolve()
adapter = root/'sweep-run_reference.py'
digest = hashlib.sha256(adapter.read_bytes()).hexdigest()
manifest = (root/'sweep-binaries.sha256').read_text().splitlines()
if not any(line.split() == [digest, 'sweep-run_reference.py'] for line in manifest):
    raise RuntimeError('N1 adapter does not match its retained SHA256 manifest')
fields = ['problem', 'criterion', 'code', 'case', 'resolution', 'threads', 'version',
          'coefficient_dofs', 'dof_convention', 'status', 'converged', 'iterations',
          'native_thread_limit', 'exit_code',
          'native_solve_wall_s', 'adapter_solve_wall_s', 'total_wall_s',
          'process_wall_s', 'peak_rss_bytes', 'B_sample_relative',
          'p_sample_relative_N1', 'flux_label_max_relative', 'iota_max_absolute',
          'force_sample_relative', 'host', 'job_id', 'compiler', 'flags',
          'adapter_sha256', 'input_sha256', 'command']
stream = sys.stdout if a.output == '-' else open(a.output, 'w', newline='')
writer = csv.DictWriter(stream, fieldnames=fields)
writer.writeheader()
for code in ['vmecpp', 'gvec']:
    for case in ['landreman-integer', 'issan-integer', 'qa']:
        for mode in ([4, 8, 12] if code == 'vmecpp' else [4, 8]):
            for threads in [1, 8]:
                out = root/'runs'/code/case/f'm{mode}-t{threads}'
                if not (out/'exit.txt').exists():
                    continue
                meta = json.loads((out/'metrics.json').read_text()) if (out/'metrics.json').exists() else {}
                errors = {}
                if (out/'oracle.txt').exists():
                    for line in (out/'oracle.txt').read_text().splitlines():
                        key, sep, value = line.partition('=')
                        if sep:
                            try:
                                errors[key] = float(value)
                            except ValueError:
                                pass
                wall, rss = '', 0
                for line in (out/'time.txt').read_text().splitlines():
                    if 'Maximum resident set size' in line:
                        rss = int(line.rsplit(':', 1)[1])*1024
                    if 'Elapsed (wall clock)' in line:
                        parts = list(map(float, line.rsplit(': ', 1)[1].split(':')))
                        wall = sum(v*60**i for i,v in enumerate(reversed(parts)))
                iota = meta.get('iota_input_convention', [])
                iota_error = max(abs(x-2) for x in iota) if iota and case != 'qa' else ''
                command = ['./mem_guard.sh', '64', '/usr/bin/time', '-v', 'timeout', '600',
                           'taskset', '-c', ','.join(map(str, meta.get('allowed_cpus', range(threads)))),
                           '/home/ert/k6bench-20261021/venv/bin/python', 'sweep-run_reference.py',
                           code, meta.get('input', 'INPUT_FROM_N1_SWEEP'), str(mode), str(out)]
                if case != 'qa':
                    command.append('--clockwise')
                writer.writerow(dict(problem='L6', criterion=4, code=code, case=case,
                    resolution=mode, threads=threads, version=meta.get('version', 'unavailable'),
                    coefficient_dofs=meta.get('coefficient_dofs', ''),
                    dof_convention=meta.get('dof_convention', ''),
                    status=meta.get('status', 'no metrics'),
                    converged=meta.get('success', meta.get('ier_flag') == 0 if code == 'vmecpp' else ''),
                    iterations=meta.get('iterations', ''), native_thread_limit=meta.get('native_thread_limit', 'unavailable'),
                    exit_code=(out/'exit.txt').read_text().strip(),
                    native_solve_wall_s=meta.get('native_solve_wall_s', ''),
                    adapter_solve_wall_s=meta.get('solve_wall_s', ''),
                    total_wall_s=meta.get('total_wall_s', ''), process_wall_s=wall,
                    peak_rss_bytes=rss, B_sample_relative=errors.get('B_sample_rms_relative', ''),
                    p_sample_relative_N1=errors.get('p_sample_rms_relative', ''),
                    flux_label_max_relative=errors.get('surface_label_max_relative', ''),
                    iota_max_absolute=iota_error,
                    force_sample_relative=errors.get('force_sample_rms_relative', ''),
                    host=meta.get('host', 'node7'), job_id=3355100,
                    compiler='external package: build compiler unavailable',
                    flags='external package: build flags unavailable', adapter_sha256=digest,
                    input_sha256=meta.get('input_sha256', ''),
                    command=f'cd {shlex.quote(str(root))} && OMP_NUM_THREADS={threads} '+shlex.join(command)))
if stream is not sys.stdout:
    stream.close()
