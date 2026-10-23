#!/usr/bin/env python3
"""Run VMEC++ on the Solov'ev case and sample it at the common L6 positions.

Usage:
  vmec_case.py --case CASEDIR --mpol 8 --ns 51 --ftol 1e-12 --threads 1 --out OUTDIR
               [--profile-type cubic_spline|akima_spline] [--boundary boundary.txt]
               [--niter 20000] [--phi-sign +1 --iota-sign +1] [--nprofile N]

CASEDIR holds profiles.txt, boundary.txt, samples.txt, reference_samples.txt and
case_info.json from solovev_reference.py.  OUTDIR receives input.json, wout.nc,
samples_out.txt (same columns as reference_samples.txt), metrics.json and l6.json.

Model: fixed boundary, nfp=1, ntor=0, lasym=False, ncurr=0 (iota-driven).  p(s) [Pa]
and iota(s) are spline tables from profiles.txt (pres_scale=1), phiedge [Wb] the exact
toroidal flux.  ftol is VMEC's force tolerance: the run stops when the invariant
force residuals fsqr, fsqz and fsql are each below ftol (checked against
wout.fsqr/fsqz/fsql/ftolv in metrics.json).

Sampling inverts (R,Z)->(s,theta) at phi=0 with Newton iterations on the Fourier
sums whose radial dependence is SIMSOPT's `vmec_splines` (cubic interpolating splines
of rmnc/zmns on the full grid); B, p, iota come from `vmec_compute_geometry` with the
same splines.  Flux psi [Wb/rad] = chi_sign*(chi(s)-chi(1))/(2 pi), chi_sign = -1,
with chi = wout.chi interpolated by the same cubic spline.  [V] VMEC++ chi is the
poloidal flux in Wb (not per radian), chi(0)=0, chi(1)=-2 pi (psi_edge-psi_axis) for
phiedge>0, iota>0, counter-clockwise boundary (wout.phi[-1] = phiedge in Wb); hence
psi(0) = chi(1)/(2 pi) = psi_axis<0 and psi(1)=0 (see metrics psi_axis_vmec).
The sign combination phi_sign=+1, iota_sign=+1 with a counter-clockwise (R,Z)
boundary reproduces the exact B_R, B_Z, B_phi, p, iota (checked at ns=25, mpol=8).
"""
import argparse
import importlib.metadata
import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np

T_START = time.perf_counter()
MU0 = 4e-7 * np.pi


def read_profiles(path):
    phiedge = None
    for line in open(path):
        if line.startswith("# phiedge ="):
            phiedge = float(line.split("=")[1])
    tab = np.loadtxt(path, comments="#")
    return phiedge, tab  # columns psiN s p iota


def read_boundary(path):
    tab = np.loadtxt(path, comments="#")
    return tab[:, 1], tab[:, 2]  # rbc_m, zbs_m for m = 0..


def build_input(case, mpol, ns, ftol, niter, profile_type, phi_sign, iota_sign,
                nprofile=None, boundary="boundary.txt"):
    import vmecpp

    phiedge, prof = read_profiles(Path(case) / "profiles.txt")
    if nprofile and nprofile < len(prof):
        idx = np.round(np.linspace(0, len(prof) - 1, nprofile)).astype(int)
        prof = prof[idx]
    rbc_m, zbs_m = read_boundary(Path(case) / boundary)
    info = json.load(open(Path(case) / "case_info.json"))
    inp = vmecpp.VmecInput(
        lasym=False, nfp=1, mpol=mpol, ntor=0,
        ntheta=0, nzeta=1,
        ns_array=np.array([ns]), ftol_array=np.array([ftol]),
        niter_array=np.array([niter]),
        phiedge=phi_sign * phiedge, ncurr=0,
        pmass_type=profile_type, am_aux_s=prof[:, 1], am_aux_f=prof[:, 2],
        pres_scale=1.0,
        piota_type=profile_type, ai_aux_s=prof[:, 1], ai_aux_f=iota_sign * prof[:, 3],
        curtor=phi_sign * iota_sign * info["i_tor_A"],
        lfreeb=False,
        raxis_c=np.array([info["r_axis"]]), zaxis_s=np.array([0.0]),
        rbc=rbc_m[:mpol].reshape(mpol, 1), zbs=zbs_m[:mpol].reshape(mpol, 1),
        return_outputs_even_if_not_converged=True,
    )
    return inp


class Sampler:
    """(R,Z)->(s,theta) at phi=0 and field evaluation on SIMSOPT's radial splines."""

    def __init__(self, wout_path):
        from simsopt.mhd import Vmec
        from simsopt.mhd.vmec_diagnostics import vmec_splines

        self.vmec = Vmec(str(wout_path))
        self.sp = vmec_splines(self.vmec)
        self.wout = self.vmec.wout
        self.xm = np.asarray(self.sp.xm, dtype=float)

    def _modes(self, s, derivative=False):
        src = (self.sp.d_rmnc_d_s, self.sp.d_zmns_d_s) if derivative else (
            self.sp.rmnc, self.sp.zmns)
        return (np.array([f(s) for f in src[0]]), np.array([f(s) for f in src[1]]))

    def rz(self, s, th):
        """R, Z, R_s, Z_s, R_th, Z_th for arrays s, th (N,)."""
        r, z = self._modes(s)
        rs, zs = self._modes(s, True)
        c, sn = np.cos(self.xm[:, None] * th[None]), np.sin(self.xm[:, None] * th[None])
        R, Z = np.sum(r * c, 0), np.sum(z * sn, 0)
        Rs, Zs = np.sum(rs * c, 0), np.sum(zs * sn, 0)
        Rt = np.sum(-r * self.xm[:, None] * sn, 0)
        Zt = np.sum(z * self.xm[:, None] * c, 0)
        return R, Z, Rs, Zs, Rt, Zt

    def invert(self, R, Z):
        sg = np.linspace(0.01, 1.0, 100)
        tg = np.linspace(0, 2 * np.pi, 360, endpoint=False)
        rm, zm = self._modes(sg)
        Rg = np.einsum("ms,mt->st", rm, np.cos(self.xm[:, None] * tg[None]))
        Zg = np.einsum("ms,mt->st", zm, np.sin(self.xm[:, None] * tg[None]))
        s = np.empty(len(R))
        th = np.empty(len(R))
        for i in range(len(R)):
            d = (Rg - R[i]) ** 2 + (Zg - Z[i]) ** 2
            a, b = np.unravel_index(np.argmin(d), d.shape)
            s[i], th[i] = sg[a], tg[b]
        for it in range(100):
            Rv, Zv, Rs, Zs, Rt, Zt = self.rz(s, th)
            fr, fz = Rv - R, Zv - Z
            det = Rs * Zt - Rt * Zs
            ds = -(fr * Zt - Rt * fz) / det
            dt = -(Rs * fz - fr * Zs) / det
            ds = np.clip(ds, -0.2, 0.2)
            s = np.clip(s + ds, 1e-6, 1.0)
            th = th + np.clip(dt, -0.5, 0.5)
            if max(np.max(np.abs(fr)), np.max(np.abs(fz))) < 1e-14:
                break
        Rv, Zv = self.rz(s, th)[:2]
        resid = float(max(np.max(np.abs(Rv - R)), np.max(np.abs(Zv - Z))))
        return s, th, resid

    def evaluate(self, R, Z, sigma_chi=1.0):
        from simsopt.mhd.vmec_diagnostics import vmec_compute_geometry

        s, th, resid = self.invert(R, Z)
        n = len(s)
        data = vmec_compute_geometry(self.sp, s, th.reshape(n, 1, 1), np.zeros((n, 1, 1)))
        g = lambda name: np.asarray(getattr(data, name)).reshape(n)  # noqa: E731
        sx, sy = g("X"), g("Y")
        assert np.max(np.abs(sx - R)) < 1e-9 and np.max(np.abs(sy)) < 1e-9
        bR, bPhi, bZ = g("B_X"), g("B_Y"), g("B_Z")
        w = self.wout
        from scipy.interpolate import InterpolatedUnivariateSpline as IUS

        sfull = np.linspace(0, 1, int(w.ns))
        chi = IUS(sfull, np.asarray(w.chi))
        psi = sigma_chi * (chi(s) - chi(1.0)) / (2 * np.pi)
        return dict(s=s, theta=th, resid=resid, R=R, Z=Z, psi=psi, B_R=bR, B_Z=bZ,
                    B_phi=bPhi, p=self.sp.pressure(s), iota=self.sp.iota(s))


def write_samples(path, d, pscale):
    absb = np.sqrt(d["B_R"] ** 2 + d["B_Z"] ** 2 + d["B_phi"] ** 2)
    with open(path, "w") as fh:
        fh.write("# R Z psi B_R B_Z B_phi absB p iota (VMEC++ via SIMSOPT; SI)\n")
        for i in range(len(d["R"])):
            fh.write(" ".join(f"{v:.17e}" for v in (
                d["R"][i], d["Z"][i], d["psi"][i], d["B_R"][i], d["B_Z"][i],
                d["B_phi"][i], absb[i], pscale * d["p"][i], d["iota"][i])) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True)
    ap.add_argument("--mpol", type=int, required=True)
    ap.add_argument("--ns", type=int, required=True)
    ap.add_argument("--ftol", type=float, required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--profile-type", default="cubic_spline")
    ap.add_argument("--boundary", default="boundary.txt")
    ap.add_argument("--niter", type=int, default=20000)
    ap.add_argument("--nprofile", type=int, default=None)
    ap.add_argument("--phi-sign", type=int, default=1)
    ap.add_argument("--iota-sign", type=int, default=1)
    ap.add_argument("--chi-sign", type=float, default=-1.0)
    ap.add_argument("--pscale", type=float, default=1.0,
                    help="multiplier on simsopt pressure (wout.pres) to give Pa")
    ap.add_argument("--no-sample", action="store_true")
    a = ap.parse_args()
    a.case = str(Path(a.case).resolve())
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(out)
    import vmecpp

    rec = dict(status="incomplete", mpol=a.mpol, ns=a.ns, ftol=a.ftol, threads=a.threads,
               profile_type=a.profile_type, boundary=a.boundary, phi_sign=a.phi_sign,
               iota_sign=a.iota_sign, version=importlib.metadata.version("vmecpp"),
               omp_num_threads=os.environ.get("OMP_NUM_THREADS"),
               allowed_cpus=sorted(os.sched_getaffinity(0)))

    def persist():
        rec["total_wall_s"] = time.perf_counter() - T_START
        rec["process_peak_rss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        (out / "metrics.json").write_text(json.dumps(rec, indent=1) + "\n")

    try:
        inp = build_input(a.case, a.mpol, a.ns, a.ftol, a.niter, a.profile_type,
                          a.phi_sign, a.iota_sign, a.nprofile, a.boundary)
        (out / "input.json").write_text(inp.model_dump_json(indent=1))
        t0 = time.perf_counter()
        res = vmecpp.run(inp, max_threads=a.threads, verbose=1)
        rec["native_solve_wall_s"] = time.perf_counter() - t0
        rec["rss_after_solve_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        w = res.wout
        w.save("wout.nc")
        mnmax = int(w.mnmax)
        rec.update(ier_flag=int(w.ier_flag), fsqr=float(w.fsqr), fsqz=float(w.fsqz),
                   fsql=float(w.fsql), ftolv=float(w.ftolv), iterations=int(w.itfsq),
                   ns_out=int(w.ns), mnmax=mnmax, coefficient_dofs=int(w.ns) * mnmax * 3,
                   native_thread_limit=min(a.threads, (int(w.ns) - 1) // 2),
                   converged=bool(int(w.ier_flag) == 0 and max(float(w.fsqr), float(w.fsqz),
                                                               float(w.fsql)) <= a.ftol),
                   phiedge_out=float(w.phi[-1]), chi_edge=float(w.chi[-1]),
                   psi_axis_vmec=float(-a.chi_sign * w.chi[-1] / (2 * np.pi)),
                   rmax_surf=float(w.rmax_surf))
        info = json.load(open(Path(a.case) / "case_info.json"))
        rec["psi_axis_exact"] = info["psi_axis"]
        rec["psi_axis_rel_err"] = rec["psi_axis_vmec"] / info["psi_axis"] - 1.0
        rec["phiedge_rel_err"] = rec["phiedge_out"] / (a.phi_sign * info["phi_edge"]) - 1.0
        persist()
        if not a.no_sample:
            smp = np.loadtxt(Path(a.case) / "samples.txt", comments="#")
            sampler = Sampler(out / "wout.nc")
            d = sampler.evaluate(smp[:, 0], smp[:, 1], a.chi_sign)
            rec["inversion_residual_m"] = d["resid"]
            rec["s_range"] = [float(d["s"].min()), float(d["s"].max())]
            write_samples(out / "samples_out.txt", d, a.pscale)
            sys.path.insert(0, str(Path(__file__).resolve().parent))
            import l6_metric

            ref = l6_metric.load(Path(a.case) / "reference_samples.txt")
            code = l6_metric.load(out / "samples_out.txt")
            span = l6_metric.psi_span(Path(a.case) / "reference_samples.txt")
            l6 = l6_metric.l6(ref, code, span)
            (out / "l6.json").write_text(json.dumps(l6) + "\n")
            rec["l6"] = l6
        rec["status"] = "ok"
    finally:
        persist()


if __name__ == "__main__":
    main()
