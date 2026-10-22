#!/usr/bin/env python3
"""Render qualified measured rows; never extrapolate accuracy targets."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import math

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('csv', nargs='+', type=Path)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
rows = []
for path in a.csv:
    rows.extend(csv.DictReader(path.open()))
rows = [r for r in rows if r['status'] == 'V' and r['dirty'] == 'False']
lines = ['# Assembled baselines', '',
         '[V] Process wall time includes assembly, constraints, preconditioner/factor setup,',
         'solve and error integration. It excludes build and Python mesh preparation.',
         'RSS is the child high-water mark from wait4/getrusage, in bytes in the CSV.',
         'Every measured row includes its exact command, source revision, executable SHA256,',
         'compiler and flags, host/job ID, mesh checksum and linear residual.', '',
         '| Problem | Criterion | Family | preconditioner | p | 1/h | Ratio | DOFs | rel L2 | axis rel L2 | chi_perp,num | iterations | wall (s) | RSS (MiB) |',
         '|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for r in rows:
    fmt = lambda k: f'{float(r[k]):.4e}' if r[k] else '—'
    lines.append(f"| {r['problem']} | {r['criterion']} | {r['family']} | {r['preconditioner']} | {r['p']} | "
                 f"{1/float(r['h']):g} | {r['anisotropy_ratio']} | {r['dofs']} | "
                 f"{fmt('rel_l2')} | {fmt('axis_rel_l2')} | {fmt('chi_perp_num')} | "
                 f"{r['iterations']} | {float(r['wall_s']):.5g} | "
                 f"{float(r['peak_rss_bytes'])/2**20:.2f} |")
lines += ['', '## Time to measured accuracy', '',
          '[V] Minimum measured process wall time among qualified rows reaching the target.',
          'A dash means no run reached it; no fitted or extrapolated timing is substituted.', '',
          '| Problem | Family | Ratio | target rel L2 | DOFs | wall (s) |',
          '|---|---|---:|---:|---:|---:|']
groups = defaultdict(list)
for r in rows:
    groups[r['problem'], r['family'], r['anisotropy_ratio']].append(r)
for key, runs in sorted(groups.items()):
    for target in (1e-3, 1e-5, 1e-7):
        hits = [r for r in runs if float(r['rel_l2']) <= target]
        hit = min(hits, key=lambda r: float(r['wall_s'])) if hits else None
        lines.append('| '+' | '.join(key)+f" | {target:.0e} | "+
                     (f"{hit['dofs']} | {float(hit['wall_s']):.5g} |" if hit else '— | — |'))
lines += ['', '## Axis refinement diagnostics', '',
          '[V] Rates use the declared 1/n mesh parameter and physical quadrature norms.',
          'These through-axis rates do not establish a rate loss relative to a polar chart.', '',
          '| Problem | Family | p | refinement | global L2 rate | axis-tube L2 rate | global minus axis |',
          '|---|---|---:|---|---:|---:|---:|']
axis_groups = defaultdict(list)
for r in rows:
    if r['axis_rel_l2'] and float(r['axis_rel_l2']) > 0:
        axis_groups[r['problem'],r['family'],r['p'],r['anisotropy_ratio'],r['preconditioner'],r['threads']].append(r)
for key, runs in sorted(axis_groups.items()):
    runs.sort(key=lambda r: -float(r['h']))
    for lo, hi in zip(runs, runs[1:]):
        if lo['h'] == hi['h']:
            continue
        denom = math.log(float(lo['h'])/float(hi['h']))
        global_rate = math.log(float(lo['rel_l2'])/float(hi['rel_l2']))/denom
        axis_rate = math.log(float(lo['axis_rel_l2'])/float(hi['axis_rel_l2']))/denom
        lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {1/float(lo['h']):g} to "
                     f"{1/float(hi['h']):g} | {global_rate:.3f} | {axis_rate:.3f} | "
                     f"{global_rate-axis_rate:.3f} |")
a.output.write_text('\n'.join(lines)+'\n')
