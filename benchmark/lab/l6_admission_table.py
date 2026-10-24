#!/usr/bin/env python3
"""Markdown tables of the L6 admission ledger written by import_n1.py --admit."""
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1])))


def num(row, key, fmt):
    return format(float(row[key]), fmt) if row[key] not in ('', None) else '—'


def decision(row):
    if row['admitted_1e-3'] == 'reference':
        return 'reference'
    if row['admitted_1e-3'] == 'True':
        return 'yes, 1e-3'
    if row['admitted_1e-2'] == 'True':
        return 'yes, 1e-2'
    return 'no' if row['e_max'] else 'no (failed)'


print('| Code | Case | m | threads | e_B | e_p | e_iota | e_flux | max | admitted |')
print('|---|---|---:|---:|---:|---:|---:|---:|---:|---|')
for r in rows:
    print(f"| {r['code']} | {r['case']} | {r['resolution']} | {r['threads']} | "
          + ' | '.join(num(r, k, '.2e') for k in
                       ('e_B', 'e_p', 'e_iota', 'e_flux', 'e_max'))
          + f' | {decision(r)} |')
print()
print('| Code | Case | m | cores | native solve (s) | adapter solve (s) '
      '| process (s) | peak RSS (MiB) | admitted |')
print('|---|---|---:|---:|---:|---:|---:|---:|---|')
for r in rows:
    if decision(r).startswith('yes'):
        print(f"| {r['code']} | {r['case']} | {r['resolution']} | {r['threads']} | "
              f"{num(r, 'native_solve_wall_s', '.4g')} | "
              f"{num(r, 'adapter_solve_wall_s', '.4g')} | "
              f"{num(r, 'process_wall_s', '.4g')} | "
              f"{float(r['peak_rss_bytes'])/2**20:.1f} | {decision(r)} |")
