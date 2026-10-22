#!/usr/bin/env python3
"""Run the authoritative KIN6D N1 reference adapter without copying its cases."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from run import ROOT, FIELDS, append_row

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--kin6d', type=Path, default=Path.home()/'code/kin6d')
p.add_argument('--code', choices=['vmecpp', 'gvec'], required=True)
p.add_argument('--resolution', type=int, required=True)
p.add_argument('--input', type=Path)
p.add_argument('--threads', type=int, required=True)
p.add_argument('--clockwise', action='store_true')
p.add_argument('--compiler', required=True)
p.add_argument('--flags', required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--csv', type=Path, default=Path('l6.csv'))
a = p.parse_args()
if a.threads < 1 or a.resolution < 1:
    p.error('positive resolution and threads required')
a.kin6d = a.kin6d.resolve()
adapter = a.kin6d/'benchmarks/nested/run_reference.py'
case = (a.input or a.kin6d/'benchmarks/nested/input.landreman-paul-qa').resolve()
command = [sys.executable, str(adapter), a.code, str(case), str(a.resolution),
           str(a.output.resolve())]+(['--clockwise'] if a.clockwise else [])
env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), OPENBLAS_NUM_THREADS='1',
           MPLBACKEND='Agg')
env.pop('DISPLAY', None)
subprocess.run(command, env=env, check=True)
metrics = json.loads((a.output/'metrics.json').read_text())
row = dict.fromkeys(FIELDS, '')
row.update(problem='L6', criterion='4', family=a.code, p='', h='',
           preconditioner='external', threads=a.threads,
           revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
           compiler=a.compiler, flags=a.flags, host=platform.node(),
           job_id=os.environ.get('SLURM_JOB_ID', 'none'), size=a.resolution,
           iterations=metrics.get('iterations', ''), wall_s=metrics['solve_wall_s'],
           peak_rss_bytes=metrics['process_peak_rss_kib']*(1 if sys.platform == 'darwin' else 1024),
           rss_method='N1.resource.getrusage.self', geometry=case.name,
           status='unqualified: matched physical accuracy not yet supplied',
           mesh_sha256=hashlib.sha256(case.read_bytes()).hexdigest(),
           command=shlex.join([sys.executable, str(Path(__file__).relative_to(ROOT))]+sys.argv[1:]))
append_row(a.csv, row)
print(json.dumps(row, indent=2))
