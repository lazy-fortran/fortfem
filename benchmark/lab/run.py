#!/usr/bin/env python3
"""Reproducible lab runs. Build once remotely, then time the binary only."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
FIELDS = ['problem', 'criterion', 'family', 'p', 'h', 'anisotropy_ratio', 'preconditioner', 'threads',
          'revision', 'executable_sha256', 'dirty', 'compiler', 'flags', 'host', 'job_id', 'size',
          'dofs', 'rel_l2', 'rel_h1_seminorm', 'rel_hcurl', 'axis_rel_l2', 'axis_ball_rel_l2', 'chi_perp_num',
          'iterations', 'wall_s', 'peak_rss_bytes', 'rss_method', 'residual',
          'nnz', 'status', 'geometry', 'mesh_sha256', 'command']


def append_row(path, row):
    path = Path(path)
    exists = path.exists() and path.stat().st_size > 0
    if exists:
        with path.open(newline='') as f:
            if next(csv.reader(f)) != FIELDS:
                raise RuntimeError('CSV schema differs; choose a new output file')
    with path.open('a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def geometry(u, v, phi):
    a = 1 + .1 * math.cos(5 * phi)
    r = 3 + .1 * math.cos(5 * phi) + a * u
    z = .1 * math.sin(5 * phi) + 1.3 * v
    return r * math.cos(phi), r * math.sin(phi), z


def determinant(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) -
            a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))


def mesh(problem, n, alignment='unaligned'):
    if problem == 'L3a':
        vertices = [(i/n-.5, j/n-.5, 0.) for j in range(n+1)
                    for i in range(n+1)]
        cells = []
        for j in range(n):
            for i in range(n):
                a = i + (n+1)*j
                cells.extend([(a, a+1, a+n+2), (a, a+n+2, a+n+1)])
        dim = 2
    else:
        # Cartesian disk triangulation: a single axis vertex, never an axis ring.
        nt, np = 8*n, 20*n
        disk = [(0., 0.)]
        for j in range(1, n+1):
            disk.extend([(j/n*math.cos(2*math.pi*i/nt),
                          j/n*math.sin(2*math.pi*i/nt)) for i in range(nt)])
        triangles = [(0, 1+i, 1+(i+1) % nt) for i in range(nt)]
        for j in range(1, n):
            for i in range(nt):
                a, b = 1+(j-1)*nt+i, 1+(j-1)*nt+(i+1) % nt
                c, d = a+nt, b+nt
                triangles.extend([(a, c, d), (a, d, b)])
        vertices = []
        for k in range(np):
            phi = 2*math.pi*k/np
            for u, v in disk:
                vertices.append(geometry(u, v, phi))
        cells = []
        nd = len(disk)
        for k in range(np):
            for tri in triangles:
                a, b, c = sorted(tri)
                a, b, c, A, B, C = (a+k*nd, b+k*nd, c+k*nd,
                                    a+(k+1) % np*nd, b+(k+1) % np*nd,
                                    c+(k+1) % np*nd)
                cells.extend([(a, b, c, C), (a, b, B, C), (a, A, B, C)])
        oriented = []
        for cell in cells:
            x = vertices[cell[0]]
            edges = [tuple(vertices[j][d]-x[d] for d in range(3)) for j in cell[1:]]
            det = determinant(*edges)
            if abs(det) < 1.e-14:
                raise ValueError('degenerate torus mesh')
            oriented.append(cell if det > 0 else (cell[0], cell[2], cell[1], cell[3]))
        cells, dim = oriented, 3
    faces = {}
    for c, cell in enumerate(cells):
        for i in range(dim+1):
            face = tuple(sorted(cell[:i]+cell[i+1:]))
            faces.setdefault(face, []).append((c, i))
    boundary = [[0]*(dim+1) for _ in cells]
    for face, uses in faces.items():
        if len(uses) == 1:
            c, i = uses[0]
            boundary[c][i] = 1
        elif len(uses) != 2:
            raise ValueError('non-manifold mesh')
    return dim, vertices, cells, boundary


def write_mesh(path, problem, n, alignment):
    dim, vertices, cells, boundary = mesh(problem, n, alignment)
    with path.open('w') as f:
        f.write(f'{dim} {len(vertices)} {len(cells)}\n')
        for point in vertices:
            f.write(' '.join(format(x, '.17g') for x in point)+'\n')
        for cell, faces in zip(cells, boundary):
            f.write(' '.join(map(str, [x+1 for x in cell]+faces))+'\n')
    return len(cells), hashlib.sha256(path.read_bytes()).hexdigest()


def run(args):
    if args.problem == 'L6':
        raise SystemExit('L6 uses KIN6D benchmarks/nested/run_reference.py; see PROBLEMS.md')
    if args.problem == 'L2' and (args.family != 'tet' or args.p != 1):
        raise SystemExit('L2 adapter requires tetrahedral Nedelec order 1')
    if args.family == 'gll-mf' and (args.problem not in ('L1', 'L3', 'L5') or
                                   args.preconditioner != 'pmg'):
        raise SystemExit('MF1 adapter supports L1/L3/L5 with pmg')
    if args.family == 'gll-polar' and (args.problem != 'L5' or args.preconditioner != 'ichol_row'):
        raise SystemExit('polar comparison requires L5 with ichol_row')
    if args.family == 'bspline' and args.problem not in ('L1', 'L5'):
        raise SystemExit('spline transport adapter pending')
    if args.alignment != 'unaligned':
        raise SystemExit('field-aligned periodic mesh closure is pending; no false aligned label')
    if args.problem == 'L3a' and args.p != 1:
        raise SystemExit('existing 2D tensor adapter supports P1')
    if args.h <= 0:
        raise SystemExit('h must be positive')
    if args.family in ('tet', 'bspline') and args.preconditioner == 'pmg':
        raise SystemExit('pmg requires the matrix-free GLL driver')
    if not args.problem.startswith('L3') and args.ratio != 1:
        raise SystemExit('anisotropy ratio applies only to L3/L3a')
    n = args.n if args.n is not None else round(1/args.h)
    if n < 1 or args.threads < 1 or args.p not in (1, 2, 3, 4):
        raise SystemExit('positive resolution/threads and degree 1..4 required')
    executable = Path(args.executable).resolve()
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads),
               OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg')
    env.pop('DISPLAY', None)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))
    compiler = subprocess.check_output([args.compiler, '--version'], text=True).splitlines()[0]
    with tempfile.TemporaryDirectory(prefix='ben1-', dir=args.scratch) as tmp:
        meshfile = Path(tmp)/'mesh.txt'
        if args.family in ('bspline', 'gll-mf', 'gll-polar'):
            size, checksum = (160 if args.family == 'gll-polar' else
                              (100 if args.family == 'gll-mf' else 80))*n**3, ''
        else:
            size, checksum = write_mesh(meshfile, args.problem, n, args.alignment)
        command = [str(executable), args.problem, str(args.p),
                   str(n) if args.family in ('bspline', 'gll-mf', 'gll-polar') else str(meshfile),
                   args.preconditioner, str(args.ratio)]
        base_row = dict(zip(FIELDS, ['']*len(FIELDS)))
        base_row.update(problem=args.problem, criterion='2' if args.problem == 'L5' else
                   ('3' if args.problem.startswith('L3') else ('4' if args.problem == 'L4' else '1')), family=args.family,
                   p=args.p, h=1/n, anisotropy_ratio=args.ratio, preconditioner=args.preconditioner, threads=args.threads,
                   revision=revision, executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(), dirty=dirty, compiler=compiler, flags=args.flags,
                   host=platform.node(), job_id=os.environ.get('SLURM_JOB_ID', 'none'),
                   size=size, dofs='', rel_l2='', rel_h1_seminorm='',
                   rel_hcurl='',
                   axis_rel_l2='',
                   axis_ball_rel_l2='', chi_perp_num='', iterations='', wall_s='',
                   peak_rss_bytes='',
                   rss_method='wait4.getrusage.child', residual='',
                   nnz='', status='pending', geometry='square' if args.problem == 'L3a'
                   else ('nfp5-polar-gll-approximation' if args.family == 'gll-polar'
                         else ('nfp5-curved-gll-approximation' if args.family == 'gll-mf'
                         else ('nfp5-periodic-spline-approximation' if args.family == 'bspline'
                         else 'nfp5-affine-tetra-approximation'))), mesh_sha256=checksum,
                   command=shlex.join([f'OMP_NUM_THREADS={args.threads}', sys.executable,
                                    str(Path(__file__).relative_to(ROOT))]+
                                   ([args.problem] if Path(sys.argv[0]).stem != 'run' else [])+sys.argv[1:]))
        records = []
        wall, peak = 0., 0
        for mode_ratio in ([3., 1.] if args.problem == 'L4' else [args.ratio]):
            command[-1] = str(mode_ratio)
            start = time.perf_counter()
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True, env=env)
            # Per-child rusage, never the maximum of preceding runs.
            pid, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
            stdout, stderr = proc.communicate()
            wall += time.perf_counter()-start
            peak = max(peak, usage.ru_maxrss)
            if proc.returncode:
                row = base_row.copy()
                row.update(status=f'failed: exit {proc.returncode}', wall_s=wall,
                           peak_rss_bytes=peak*(1 if sys.platform == 'darwin' else 1024))
                append_row(args.output, row)
                sys.stderr.write(''.join((stdout+stderr).splitlines(True)[-100:]))
                raise SystemExit(proc.returncode)
            record = [line.split()[1:] for line in stdout.splitlines()
                      if line.startswith('LAB ')]
            if len(record) != 1:
                raise RuntimeError('expected exactly one LAB record per modal solve')
            records.append(record[0])
        dofs, l2, h1, axis, center, residual, nnz, iterations = records[0][:8]
        ball = float(records[0][8]) if len(records[0]) == 9 else ''
        if args.problem == 'L4':
            # Constant SPD block diagonalises into eigenvalues 3,1,1.
            # Exact component amplitudes in this orthonormal frame have
            # squared magnitudes 4.5,0.5,9; the eigenvalue-1 solve is reused.
            other = records[1]
            dofs = 3*int(dofs)
            l2 = math.sqrt((4.5*float(l2)**2+9.5*float(other[1])**2)/14)
            h1 = math.sqrt((4.5*float(h1)**2+9.5*float(other[2])**2)/14)
            axis = math.sqrt((4.5*float(axis)**2+9.5*float(other[3])**2)/14)
            residual = max(float(residual), float(other[5]))
            nnz = int(float(nnz))+int(float(other[6]))
            iterations = int(iterations)+int(other[7])
        pollution = ''
        if args.problem in ('L3a', 'L3'):
            if float(center) <= 0:
                raise RuntimeError('nonpositive central temperature')
            # chi_perp=1; cosine test has T_exact(0,0)=1.
            pollution = 1/float(center)-1
        row = base_row.copy()
        row.update(dofs=int(dofs), rel_l2=float(l2),
                   rel_h1_seminorm=float(h1) if args.problem != 'L2' else '',
                   rel_hcurl=float(h1) if args.problem == 'L2' else '',
                   axis_rel_l2=float(axis) if args.problem not in ('L3a', 'L2') else '',
                   axis_ball_rel_l2=ball, chi_perp_num=pollution,
                   iterations=int(iterations), wall_s=wall,
                   peak_rss_bytes=peak*(1 if sys.platform == 'darwin' else 1024),
                   residual=float(residual), nnz=int(float(nnz)), status='V')
        append_row(args.output, row)
        print(json.dumps(row, indent=2))


def parser(problem=None):
    p = argparse.ArgumentParser(description=__doc__)
    if problem is None:
        p.add_argument('problem', choices=['L1', 'L2', 'L3', 'L3a', 'L4', 'L5', 'L6'])
    else:
        p.set_defaults(problem=problem)
    p.add_argument('--family', choices=['tet', 'bspline', 'gll-mf', 'gll-polar'], default='tet')
    p.add_argument('--p', type=int, default=1)
    resolution = p.add_mutually_exclusive_group()
    resolution.add_argument('--n', type=int)
    resolution.add_argument('--h', type=float, default=1.)
    p.add_argument('--preconditioner', choices=['ldlt', 'ichol', 'ichol_row', 'ilu', 'jacobi', 'pmg'], default='ichol_row')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--ratio', type=float, choices=[1., 100., 10000., 1000000.], default=1.)
    p.add_argument('--alignment', choices=['unaligned', 'aligned'], default='unaligned')
    p.add_argument('--executable', required=True)
    p.add_argument('--compiler', default='gfortran')
    p.add_argument('--flags', required=True, help='exact flags used to build executable')
    p.add_argument('--scratch', default=os.environ.get('TMPDIR', '/tmp'))
    p.add_argument('--output', default='lab.csv')
    return p


if __name__ == '__main__':
    run(parser().parse_args())
