#!/usr/bin/env python3
"""Read the existing N1 sweep ledger; never rerun or copy equilibrium cases.

With --admit the L6 matched-accuracy vector (e_B, e_p, e_iota, e_flux) of
PROBLEMS.md is evaluated for every successful record against the finest VMEC++
reference of its case, by reading the retained samples.txt and wout files only.
That mode needs numpy, scipy and SIMSOPT (the N1 virtual environment)."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shlex
import sys

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('root', type=Path)
p.add_argument('--output', required=True)
p.add_argument('--admit', action='store_true',
               help='add the L6 admission columns (needs the N1 SIMSOPT environment)')
a = p.parse_args()
root = a.root.resolve()
adapter = root/'sweep-run_reference.py'
digest = hashlib.sha256(adapter.read_bytes()).hexdigest()
manifest = (root/'sweep-binaries.sha256').read_text().splitlines()
if not any(line.split() == [digest, 'sweep-run_reference.py'] for line in manifest):
    raise RuntimeError('N1 adapter does not match its retained SHA256 manifest')
# Finest VMEC++ run per case. QA uses N1's m16/n12/ns201 reference, which is not
# one of the thirty ledger records. The integer families use their m12 record.
REFERENCES = {'qa': 'vmecpp/qa/m16-n12-ns201-original-t8',
              'landreman-integer': 'vmecpp/landreman-integer/m12-t1',
              'issan-integer': 'vmecpp/issan-integer/m12-t1'}
ADMISSION_FIELDS = ['reference_run', 'reference_wout_sha256', 'admission_samples',
                    'position_inverse_max_m', 'e_B', 'e_p', 'e_iota', 'e_flux', 'e_max',
                    'admitted_1e-2', 'admitted_1e-3', 'admission_note']


def invert_samples(vs, mu0, rows):
    """Reference (s, theta) at the physical sample positions, then B, p and iota there.

    Same batched inverse as N1's compare_qa_final.py: R and Z are inverted at the
    physical toroidal angle with a sparse Gauss-Newton solve, and the forward
    positions are checked afterwards."""
    import numpy as np
    from scipy.optimize import least_squares
    from scipy.sparse import coo_matrix
    from simsopt.mhd.vmec_diagnostics import vmec_compute_geometry

    def geometry(s, theta, phi):
        return vmec_compute_geometry(vs, np.atleast_1d(s),
                                     np.asarray(theta).reshape(-1, 1, 1),
                                     np.asarray(phi).reshape(-1, 1, 1))

    def vector(data, *keys):
        return np.column_stack([np.asarray(getattr(data, k)).ravel() for k in keys])

    xyz = rows[:, 1:4]
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    radius = np.hypot(xyz[:, 0], xyz[:, 1])
    theta_grid = np.arange(32)*2*np.pi/32
    count = len(rows)
    start = np.zeros((count, 2))
    guesses = {}
    for k in range(count):
        key = (float(rows[k, 0]), round(float(phi[k]), 12))
        if key not in guesses:
            guesses[key] = vector(vmec_compute_geometry(vs, [rows[k, 0]], theta_grid,
                                                        [phi[k]]), 'X', 'Y', 'Z')
        near = np.argmin(np.linalg.norm(guesses[key] - xyz[k], axis=1))
        start[k] = [rows[k, 0], theta_grid[near]]
    index = np.arange(count)
    jr = np.column_stack((2*index, 2*index, 2*index + 1, 2*index + 1)).ravel()
    jc = np.column_stack((2*index, 2*index + 1, 2*index, 2*index + 1)).ravel()
    cache = {}

    def evaluate(flat):
        if 'x' not in cache or not np.array_equal(flat, cache['x']):
            st = flat.reshape(-1, 2)
            q = geometry(st[:, 0], st[:, 1], phi)
            cache['x'] = flat.copy()
            cache['f'] = np.column_stack((q.R.ravel() - radius,
                                          q.Z.ravel() - xyz[:, 2])).ravel()
            data = np.column_stack((q.d_R_d_s.ravel(), q.d_R_d_theta_vmec.ravel(),
                                    q.d_Z_d_s.ravel(),
                                    q.d_Z_d_theta_vmec.ravel())).ravel()
            cache['j'] = coo_matrix((data, (jr, jc)), shape=(2*count, 2*count)).tocsr()
        return cache

    solution = least_squares(lambda x: evaluate(x)['f'], start.ravel(),
                             jac=lambda x: evaluate(x)['j'],
                             bounds=(np.tile([1e-10, -np.inf], count),
                                     np.tile([1 - 1e-10, np.inf], count)),
                             ftol=1e-13, xtol=1e-13, gtol=1e-13, max_nfev=50,
                             tr_options=dict(atol=1e-14, btol=1e-14))
    coordinates = solution.x.reshape(-1, 2)
    g = geometry(coordinates[:, 0], coordinates[:, 1], phi)
    position_error = float(np.max(np.linalg.norm(vector(g, 'X', 'Y', 'Z') - xyz,
                                                 axis=1)))
    if position_error > 1e-9:
        raise RuntimeError(f'forward position check failed: {position_error}')
    return (coordinates[:, 0], vector(g, 'B_X', 'B_Y', 'B_Z'),
            mu0*vs.pressure(coordinates[:, 0]), position_error)


class AdmissionReference:
    """One finest-VMEC++ reference: its wout splines and input-convention iota."""

    def __init__(self, root, name):
        from simsopt.mhd import Vmec
        from simsopt.mhd.vmec_diagnostics import vmec_splines
        self.name = name
        self.wout = root/'runs'/name/'wout_reference.nc'
        self.wout_sha256 = hashlib.sha256(self.wout.read_bytes()).hexdigest()
        self.splines = vmec_splines(Vmec(str(self.wout)))
        metrics = json.loads((root/'runs'/name/'metrics.json').read_text())
        self.iota = metrics['iota_input_convention']
        self.mu0 = 4e-7*math.pi

    def vector(self, record_dir, converged):
        """(e_B, e_p, e_iota, e_flux) of one record; equal sample weights."""
        import numpy as np
        rows = np.loadtxt(record_dir/'samples.txt')
        if rows.shape != (1728, 14) or not np.all(np.isfinite(rows)):
            raise ValueError(f'invalid physical samples: {record_dir}')
        label, b, p_ref, position_error = invert_samples(self.splines, self.mu0, rows)
        meta = json.loads((record_dir/'metrics.json').read_text())
        iota = np.asarray(meta['iota_input_convention'])
        iota_ref = np.asarray(self.iota)
        def rms(x):
            return float(np.sqrt(np.mean(np.square(x))))

        b_ref_rms = rms(np.linalg.norm(b, axis=1))
        e = dict(position_inverse_max_m=position_error,
                 e_B=rms(rows[:, 4:7] - b)/rms(b),
                 e_p=rms(rows[:, 7] - p_ref)/max(rms(p_ref), b_ref_rms**2/2),
                 e_iota=rms(iota - iota_ref)/max(rms(iota_ref), 1),
                 e_flux=rms(rows[:, 0] - label))
        e['e_max'] = max(e[k] for k in ('e_B', 'e_p', 'e_iota', 'e_flux'))
        e['admitted_1e-2'] = bool(converged and e['e_max'] <= 1e-2)
        e['admitted_1e-3'] = bool(converged and e['e_max'] <= 1e-3)
        return e

fields = ['problem', 'criterion', 'code', 'case', 'resolution', 'threads', 'version',
          'coefficient_dofs', 'dof_convention', 'status', 'converged', 'iterations',
          'native_thread_limit', 'exit_code',
          'native_solve_wall_s', 'adapter_solve_wall_s', 'total_wall_s',
          'process_wall_s', 'peak_rss_bytes', 'B_sample_relative',
          'p_sample_relative_N1', 'flux_label_max_relative', 'iota_max_absolute',
          'force_sample_relative', 'host', 'job_id', 'compiler', 'flags',
          'adapter_sha256', 'input_sha256', 'command']
if a.admit:
    fields = fields[:-1] + ADMISSION_FIELDS + ['command']
references = {}


def admission_columns(root, out, case, converged):
    """Admission vector of one record, or a note when it has no samples."""
    if not (out/'samples.txt').exists() or (out/'exit.txt').read_text().strip() != '0':
        return dict(admission_note='failed record: no samples, not admitted')
    name = REFERENCES[case]
    if name not in references:
        references[name] = AdmissionReference(root, name)
    ref = references[name]
    e = ref.vector(out, converged)
    e.update(reference_run=name, reference_wout_sha256=ref.wout_sha256,
             admission_samples=1728,
             admission_note='finest VMEC++ reference itself' if
             out == root/'runs'/name else '')
    return e


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
                row = dict(problem='L6', criterion=4, code=code, case=case,
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
                    command=f'cd {shlex.quote(str(root))} && OMP_NUM_THREADS={threads} '+shlex.join(command))
                converged = row['converged'] is True or str(row['converged']) == 'True'
                row.update(admission_columns(root, out, case, converged) if a.admit else {})
                writer.writerow(row)
if stream is not sys.stdout:
    stream.close()
