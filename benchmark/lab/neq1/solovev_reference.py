#!/usr/bin/env python3
"""Exact Solov'ev reference (Cerfon & Freidberg 2010, Phys. Plasmas 17, 032502).

Normalised x=R/R0, y=Z/R0, psi = psi0*pb(x,y),
pb = x^4/8 + A (x^2 ln x/2 - x^4/8) + sum_{i=1..7} c_i psi_i,
Delta* pb = (1-A) x^2 + A, plasma = {pb<0} bounded by pb=0 around the axis.
Physical scaling R0 [m], F0 [T m], mu0 = 4e-7 pi.  psi = poloidal flux per radian
[Wb/rad], psi_edge = 0, psi_axis < 0.  Toroidal flux Phi(psi) = 2 pi int q dpsi [Wb]
(Phi = int B_phi dR dZ, no further 2 pi).  Writes case.nml, samples.txt,
reference_samples.txt, profiles.txt, boundary.txt (see main()).

Usage: solovev_reference.py OUTDIR [--q0 1.2] [--eps .32 --kappa 1.7 --delta .33 --A -.155]
"""
import argparse
import math
import os

import numpy as np
import sympy as sp
from numpy.polynomial import chebyshev as C
from scipy.optimize import brentq

MU0 = 4e-7 * math.pi
X, Y = sp.symbols("x y", positive=True)
LN = sp.log(X)
HOMOGENEOUS = [
    sp.Integer(1),
    X**2,
    Y**2 - X**2 * LN,
    X**4 - 4 * X**2 * Y**2,
    2 * Y**4 - 9 * Y**2 * X**2 + 3 * X**4 * LN - 12 * X**2 * Y**2 * LN,
    X**6 - 12 * X**4 * Y**2 + 8 * X**2 * Y**4,
    8 * Y**6 - 140 * Y**4 * X**2 + 75 * Y**2 * X**4 - 15 * X**6 * LN
    + 180 * X**4 * Y**2 * LN - 120 * X**2 * Y**4 * LN,
]


def gs_operator(f):
    """Delta* f = x d/dx(f_x/x) + f_yy."""
    return X * sp.diff(sp.diff(f, X) / X, X) + sp.diff(f, Y, 2)


def particular(a):
    return X**4 / 8 + a * (X**2 * LN / 2 - X**4 / 8)


class Solovev:
    def __init__(self, eps=0.32, kappa=1.7, delta=0.33, A=-0.155, f0=1.0, r0=1.0,
                 q0=1.2):
        self.eps, self.kappa, self.delta, self.A = eps, kappa, delta, A
        self.f0, self.r0 = f0, r0
        # symbolic derivation: homogeneous residuals vanish, particular solves GS
        for h in HOMOGENEOUS:
            assert sp.simplify(gs_operator(h)) == 0
        ar = sp.nsimplify(A)  # exact rational for the symbolic identity check
        assert sp.simplify(gs_operator(particular(ar)) - ((1 - ar) * X**2 + ar)) == 0
        self._fun = {}
        self._build_basis()
        self.c = self._solve_bc()
        pb = particular(A) + sum(ci * h for ci, h in zip(self.c_sym, HOMOGENEOUS))
        self._make_numeric(pb)
        self._find_axis()
        self.q0 = q0
        self._set_psi0(q0)

    # ---- symbolic -> numeric
    def _build_basis(self):
        funcs = [particular(self.A)] + HOMOGENEOUS
        self.basis_syms = {
            name: [sp.diff(f, *args) if args else f for f in funcs]
            for name, args in dict(v=(), x=(X,), y=(Y,), xx=(X, X), yy=(Y, Y),
                                   xy=(X, Y)).items()}
        self.basis_num = {k: sp.lambdify((X, Y), v, "numpy")
                          for k, v in self.basis_syms.items()}

    def _row(self, name, x, y):
        return np.array([float(g) for g in self.basis_num[name](x, y)])

    def _solve_bc(self):
        e, k, d = self.eps, self.kappa, self.delta
        al = math.asin(d)
        n1 = -(1 + al) ** 2 / (e * k**2)
        n2 = (1 - al) ** 2 / (e * k**2)
        n3 = -k / (e * math.cos(al) ** 2)
        xo, xi, xt, yt = 1 + e, 1 - e, 1 - d * e, k * e
        rows, rhs = [], []

        def add(vec, particular_val):
            rows.append(vec[1:])
            rhs.append(-particular_val)

        r = lambda n, x, y: self._row(n, x, y)  # noqa: E731
        for (x, y) in ((xo, 0.0), (xi, 0.0), (xt, yt)):
            v = r("v", x, y)
            add(v, v[0])
        v = r("x", xt, yt)
        add(v, v[0])
        v = r("yy", xo, 0.0) + n1 * r("x", xo, 0.0)
        add(v, v[0])
        v = r("yy", xi, 0.0) + n2 * r("x", xi, 0.0)
        add(v, v[0])
        v = r("xx", xt, yt) + n3 * r("y", xt, yt)
        add(v, v[0])
        mat = np.array(rows)
        self.bc_cond = np.linalg.cond(mat)
        c = np.linalg.solve(mat, np.array(rhs))
        self.c_sym = list(c)
        self.bc_residual = np.max(np.abs(mat @ c - np.array(rhs)))
        return c

    def _make_numeric(self, pb):
        self.f = {k: sp.lambdify((X, Y), sp.diff(pb, *a) if a else pb, "numpy")
                  for k, a in dict(v=(), x=(X,), y=(Y,), xx=(X, X), yy=(Y, Y),
                                   xy=(X, Y)).items()}

    def _find_axis(self):
        self.xa = brentq(lambda x: self.f["x"](x, 0.0), 0.7, 1.3, xtol=1e-15, rtol=1e-15)
        self.ya = 0.0
        self.pb_a = float(self.f["v"](self.xa, 0.0))
        self.pbxx = float(self.f["xx"](self.xa, 0.0))
        self.pbyy = float(self.f["yy"](self.xa, 0.0))
        assert self.pb_a < 0 and self.pbxx > 0 and self.pbyy > 0
        assert abs(self.f["y"](self.xa, 0.0)) < 1e-14

    def _set_psi0(self, q0):
        a, k2 = self.A, self.pbxx * self.pbyy
        den = q0**2 * self.xa**2 * k2 + 2 * a * self.pb_a
        assert den > 0
        # q_axis = R0 F_a / (x_a psi0 K), F_a^2 = F0^2 - 2 A psi0^2 pb_a / R0^2
        self.psi0 = math.sqrt(self.f0**2 / den) * self.r0
        self.psi_axis = self.psi0 * self.pb_a
        self.p_axis = -self.psi0**2 * (1 - a) * self.pb_a / (MU0 * self.r0**4)
        self.cf = -2 * a * self.psi0 * self.psi_axis / self.r0**2

    # ---- exact fields (R,Z in metres)
    def psibar(self, R, Z):
        return self.f["v"](R / self.r0, Z / self.r0)

    def psi(self, R, Z):
        return self.psi0 * self.psibar(R, Z)

    def psiN(self, R, Z):
        return 1.0 - self.psibar(R, Z) / self.pb_a

    def F(self, R, Z):
        return np.sqrt(self.f0**2 - 2 * self.A * self.psi0**2 * self.psibar(R, Z)
                       / self.r0**2)

    def pressure(self, R, Z):
        return -self.psi0**2 * (1 - self.A) / (MU0 * self.r0**4) * self.psibar(R, Z)

    def field(self, R, Z):
        x, y = R / self.r0, Z / self.r0
        dpsi_dR = self.psi0 * self.f["x"](x, y) / self.r0
        dpsi_dZ = self.psi0 * self.f["y"](x, y) / self.r0
        return -dpsi_dZ / R, dpsi_dR / R, self.F(R, Z) / R

    # ---- ray geometry about the axis (polar angle t)
    def ray_root(self, c, t):
        """r(t) with pb(axis + r e(t)) = c, vectorised hybrid Newton/bisection.

        c is normalised (pb units), in [pb_a, 0], scalar or array shaped like t.
        """
        t = np.asarray(t, dtype=float)
        c = np.broadcast_to(np.asarray(c, dtype=float), t.shape)
        ct, st = np.cos(t), np.sin(t)
        lo = np.zeros_like(t)
        hi = np.full_like(t, 0.7)
        fhi = self.f["v"](self.xa + hi * ct, self.ya + hi * st) - c
        assert np.all(fhi > 0), "bracket failure (ray leaves plasma)"
        r = np.sqrt(np.clip((c - self.pb_a) / (-self.pb_a), 0, None)) * 0.3 + 1e-3
        r = np.minimum(r, 0.69)
        for _ in range(200):
            xx, yy = self.xa + r * ct, self.ya + r * st
            fv = self.f["v"](xx, yy) - c
            lo = np.where(fv < 0, r, lo)
            hi = np.where(fv >= 0, r, hi)
            g = self.f["x"](xx, yy) * ct + self.f["y"](xx, yy) * st
            with np.errstate(divide="ignore", invalid="ignore"):
                rn = r - fv / g
            bad = ~np.isfinite(rn) | (rn <= lo) | (rn >= hi)
            rn = np.where(bad, 0.5 * (lo + hi), rn)
            if np.max(np.abs(rn - r)) < 1e-16:
                r = rn
                break
            r = rn
        return r

    def contour(self, psin, t):
        """(R, Z) on the surface psiN=psin at polar angles t (metres)."""
        c = self.pb_a * (1.0 - psin)
        r = self.ray_root(c, t)
        return self.r0 * (self.xa + r * np.cos(t)), self.r0 * (self.ya + r * np.sin(t))

    def boundary_miller(self, t):
        """Edge contour on the Miller rays x=1+rho*eps*cos(t+al sin t), y=rho*kappa*eps*sin t.

        Same shape as the Miller-type parametrisation used for the CF boundary
        conditions, hence spectrally much faster in m than the polar angle.
        Returns (R, Z) in metres.
        """
        t = np.asarray(t, dtype=float)
        al = math.asin(self.delta)
        ux = self.eps * np.cos(t + al * np.sin(t))
        uy = self.kappa * self.eps * np.sin(t)
        lo, hi = np.full_like(t, 0.3), np.full_like(t, 1.8)
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            inside = self.f["v"](1.0 + mid * ux, mid * uy) < 0
            lo = np.where(inside, mid, lo)
            hi = np.where(inside, hi, mid)
        rho = 0.5 * (lo + hi)
        return self.r0 * (1.0 + rho * ux), self.r0 * rho * uy

    def q_of_psin(self, psin, nt=256):
        """Safety factor of surface psiN in (0,1] (ray quadrature, spectral)."""
        t = 2 * np.pi * np.arange(nt) / nt
        c = self.pb_a * (1.0 - psin)
        r = self.ray_root(c, t)
        x, y = self.xa + r * np.cos(t), self.ya + r * np.sin(t)
        # normalised: dl/(R|grad psi|) = r dt /(x * grad_pb.rhat) * R0^2 / (psi0 R0)
        g = self.f["x"](x, y) * np.cos(t) + self.f["y"](x, y) * np.sin(t)
        integral = np.mean(r / (x * g)) * 2 * np.pi          # normalised
        # physical: dl = R0 dl_n, R = R0 x, |grad psi| = psi0 |grad_n pb|/R0
        # dl/(R|grad psi|) = R0 * R0 /(R0 * psi0) * (dl_n/(x |grad pb|)) = R0/psi0 * ...
        Fs = math.sqrt(self.f0**2 - 2 * self.A * self.psi0**2 * c / self.r0**2)
        return Fs / (2 * np.pi) * integral * self.r0 / self.psi0

    def q_axis(self):
        return self.f0 * math.sqrt(1 - 2 * self.A * self.psi0**2 * self.pb_a
                                   / (self.f0**2 * self.r0**2)) * self.r0 / (
            self.xa * self.psi0 * math.sqrt(self.pbxx * self.pbyy))

    def build_profiles(self, n=64):
        """Chebyshev-Lobatto fit of q(psiN); Phi(psiN) via 2 pi int q dpsi."""
        j = np.arange(n + 1)
        nodes = 0.5 * (1 - np.cos(np.pi * j / n))             # psiN in [0,1]
        q = np.empty(n + 1)
        q[0] = self.q_axis()
        for i in range(1, n + 1):
            q[i] = self.q_of_psin(nodes[i])
        xs = 2 * nodes - 1
        self.q_cheb = C.chebfit(xs, q, n)
        # Phi(psiN) = 2 pi * (psi_edge - psi_axis) * int_0^psiN q dpsiN
        dpsi = -self.psi_axis
        integ = C.chebint(self.q_cheb, lbnd=-1.0)             # in xs, d xs = 2 dpsiN
        self.phi_cheb = 2 * np.pi * dpsi * 0.5 * integ
        self.q_nodes = (nodes, q)
        self.phi_edge = float(C.chebval(1.0, self.phi_cheb))

    def q_fit(self, psin):
        return C.chebval(2 * np.asarray(psin) - 1, self.q_cheb)

    def phi_fit(self, psin):
        return C.chebval(2 * np.asarray(psin) - 1, self.phi_cheb)

    def s_of_psin(self, psin):
        return self.phi_fit(psin) / self.phi_edge

    def psin_of_s(self, s):
        s = np.atleast_1d(s)
        return np.array([0.0 if v <= 0 else 1.0 if v >= 1 else
                         brentq(lambda p: self.s_of_psin(p) - v, 0.0, 1.0,
                                xtol=1e-16, rtol=1e-15) for v in s])

    def phi_direct(self, nt=256, nr=48):
        """Independent Phi_edge = int F/R dR dZ (polar Gauss-Legendre)."""
        t = 2 * np.pi * np.arange(nt) / nt
        re = self.ray_root(0.0, t)
        u, w = np.polynomial.legendre.leggauss(nr)
        tot = 0.0
        for ui, wi in zip(u, w):
            rho = 0.5 * (ui + 1) * re
            x, y = self.xa + rho * np.cos(t), self.ya + rho * np.sin(t)
            tot += wi * 0.5 * np.sum(
                self.F(self.r0 * x, self.r0 * y) / (self.r0 * x) * rho * re) * self.r0**2
        return tot * 2 * np.pi / nt

    def volume_average_p(self, nt=256, nr=48):
        t = 2 * np.pi * np.arange(nt) / nt
        re = self.ray_root(0.0, t)
        u, w = np.polynomial.legendre.leggauss(nr)
        pv = vol = 0.0
        for ui, wi in zip(u, w):
            rho = 0.5 * (ui + 1) * re
            x, y = self.xa + rho * np.cos(t), self.ya + rho * np.sin(t)
            dv = wi * 0.5 * np.sum(x * rho * re)
            vol += dv
            pv += wi * 0.5 * np.sum(self.pressure(self.r0 * x, self.r0 * y) * x
                                    * rho * re)
        return pv / vol

    def toroidal_current(self, nt=2048):
        """I_tor = (1/mu0) oint B_p . dl (counter-clockwise (R,Z) contour) [A]."""
        t = 2 * np.pi * np.arange(nt) / nt
        r = self.ray_root(0.0, t)
        # d(R,Z)/dt via analytic derivative of r(t): use spectral differentiation
        k = np.fft.fftfreq(nt, 1.0 / nt)
        dr = np.real(np.fft.ifft(1j * k * np.fft.fft(r)))
        dR = (dr * np.cos(t) - r * np.sin(t)) * self.r0
        dZ = (dr * np.sin(t) + r * np.cos(t)) * self.r0
        R, Z = self.r0 * (self.xa + r * np.cos(t)), self.r0 * r * np.sin(t)
        bR, bZ, _ = self.field(R, Z)
        return np.mean(bR * dR + bZ * dZ) * 2 * np.pi / MU0

    def beta_tor(self):
        b0 = self.f0 / self.r0
        return 2 * MU0 * self.volume_average_p() / b0**2


def fortran_real(v):
    return repr(float(v))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outdir")
    ap.add_argument("--q0", type=float, default=1.2)
    ap.add_argument("--eps", type=float, default=0.32)
    ap.add_argument("--kappa", type=float, default=1.7)
    ap.add_argument("--delta", type=float, default=0.33)
    ap.add_argument("--A", type=float, default=-0.155)
    ap.add_argument("--nprofile", type=int, default=801,
                    help="rows of profiles.txt (uniform in s)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda n: os.path.join(args.outdir, n)  # noqa: E731
    S = Solovev(args.eps, args.kappa, args.delta, args.A, q0=args.q0)
    S.build_profiles()

    # --- case.nml
    with open(out("case.nml"), "w") as fh:
        fh.write("&solovev\n")
        for k, v in (("r0", S.r0), ("f0", S.f0), ("psi0", S.psi0), ("mu0", MU0),
                     ("a", S.A), ("eps", S.eps), ("kappa", S.kappa),
                     ("delta", S.delta)):
            fh.write(f"  {k} = {fortran_real(v)}\n")
        fh.write("  c = " + ", ".join(fortran_real(v) for v in S.c) + "\n")
        for k, v in (("psi_axis", S.psi_axis), ("r_axis", S.r0 * S.xa),
                     ("z_axis", 0.0), ("p_axis", S.p_axis), ("cf", S.cf),
                     ("alpha_p", 1.0), ("alpha_f", 1.0)):
            fh.write(f"  {k} = {fortran_real(v)}\n")
        fh.write("/\n")

    # --- samples
    s_lab = [0.1, 0.25, 0.5, 0.75, 0.9]
    ang = 2 * np.pi * (np.arange(24) + 0.37) / 24
    pts = []
    psin_s = S.psin_of_s(np.array(s_lab))
    for sv, pn in zip(s_lab, psin_s):
        R, Z = S.contour(pn, ang)
        for a_, r_, z_ in zip(ang, R, Z):
            pts.append((sv, pn, r_, z_))
    with open(out("samples.txt"), "w") as fh:
        fh.write("# Solov'ev exact flux-surface sample positions: R Z [m] at phi=0\n")
        fh.write("# s = Phi/Phi_edge in {0.1,0.25,0.5,0.75,0.9}; 24 polar angles about the\n")
        fh.write("# magnetic axis 2*pi*(k+0.37)/24; point order: s-major, angle-minor\n")
        for sv, pn, r_, z_ in pts:
            fh.write(f"{r_:.17e} {z_:.17e}\n")
    P = np.array([(p[2], p[3]) for p in pts])
    sv_arr = np.array([p[0] for p in pts])
    psin = S.psiN(P[:, 0], P[:, 1])
    bR, bZ, bP = S.field(P[:, 0], P[:, 1])
    absb = np.sqrt(bR**2 + bZ**2 + bP**2)
    iota = 1.0 / S.q_fit(psin)
    with open(out("reference_samples.txt"), "w") as fh:
        fh.write("# R Z psi B_R B_Z B_phi absB p iota (exact Solov'ev; SI: m, Wb/rad, T, Pa)\n")
        fh.write(f"# psi_axis = {S.psi_axis!r}  psi_edge = 0.0  phi_edge = {S.phi_edge!r}\n")
        for i in range(len(P)):
            fh.write(" ".join(f"{v:.17e}" for v in (
                P[i, 0], P[i, 1], S.psi(P[i, 0], P[i, 1]), bR[i], bZ[i], bP[i],
                absb[i], S.pressure(P[i, 0], P[i, 1]), iota[i])) + "\n")

    # --- profiles uniform in s
    s_tab = np.linspace(0.0, 1.0, args.nprofile)
    pn_tab = S.psin_of_s(s_tab)
    p_tab = S.p_axis * (1 - pn_tab)
    iota_tab = 1.0 / S.q_fit(pn_tab)
    with open(out("profiles.txt"), "w") as fh:
        fh.write(f"# phiedge = {S.phi_edge!r}\n")
        fh.write("# phiedge: toroidal flux through the plasma cross-section [Wb] =\n")
        fh.write("#   int B_phi dR dZ = 2*pi*int q dpsi, psi = poloidal flux per radian [Wb/rad]\n")
        fh.write("# columns: psiN s p[Pa] iota ; rows uniform in s=Phi/phiedge\n")
        for a, b, c, d in zip(pn_tab, s_tab, p_tab, iota_tab):
            fh.write(f"{a:.17e} {b:.17e} {c:.17e} {d:.17e}\n")

    # --- boundary (Miller-ray parametrisation = boundary.txt; polar = boundary_polar.txt)
    nt = 1024
    t = 2 * np.pi * np.arange(nt) / nt
    td = 2 * np.pi * (np.arange(4096) + 0.5) / 4096
    mmax = 32
    kinds = {"boundary.txt": ("miller", S.boundary_miller),
             "boundary_polar.txt": ("polar", lambda tt: S.contour(1.0, tt))}
    for fname, (kind, fun) in kinds.items():
        Rb, Zb = fun(t)
        rbc = np.array([(1.0 if m == 0 else 2.0) * np.mean(Rb * np.cos(m * t))
                        for m in range(mmax + 1)])
        zbs = np.array([2.0 * np.mean(Zb * np.sin(m * t)) for m in range(mmax + 1)])
        Rd, Zd = fun(td)
        with open(out(fname), "w") as fh:
            if kind == "miller":
                fh.write("# Boundary psibar=0 contour, Miller-ray parametrisation t:\n")
                fh.write("# (R,Z)=(1+rho(t) eps cos(t+asin(delta) sin t), rho(t) kappa eps sin t)*R0\n")
            else:
                fh.write("# Boundary psibar=0 contour, polar angle t about the magnetic axis:\n")
            fh.write("# R(t) = sum_m rbc_m cos(m t), Z(t) = sum_m zbs_m sin(m t) [m];\n")
            fh.write("# VMEC theta=t, counter-clockwise in the (R,Z) plane.\n")
            fh.write("# max reconstruction error (m) max_t |(R,Z)_trunc - (R,Z)_exact| for m<mpol:\n")
            for mp in (4, 6, 8, 10, 12, 16, 24):
                ms = np.arange(mp)
                Rt = sum(rbc[m] * np.cos(m * td) for m in ms)
                Zt = sum(zbs[m] * np.sin(m * td) for m in ms)
                err = float(np.max(np.hypot(Rt - Rd, Zt - Zd)))
                fh.write(f"# trunc mpol={mp} max_err={err:.6e}\n")
            fh.write("# columns: m rbc zbs\n")
            for m in range(mmax + 1):
                fh.write(f"{m} {rbc[m]:.17e} {zbs[m]:.17e}\n")
    Rd, Zd = S.contour(1.0, td)

    # --- diagnostics to stdout (and verification of internal consistency)
    xs = np.linspace(0, 1, 5)
    ext = dict(Rmax=Rd.max(), Rmin=Rd.min(), Zmax=Zd.max())
    eps_n = (ext["Rmax"] - ext["Rmin"]) / (ext["Rmax"] + ext["Rmin"])
    kap_n = ext["Zmax"] / (0.5 * (ext["Rmax"] - ext["Rmin"]))
    iz = np.argmax(Zd)
    del_n = (0.5 * (ext["Rmax"] + ext["Rmin"]) - Rd[iz]) / (0.5 * (ext["Rmax"] - ext["Rmin"]))
    print("c =", [float(v) for v in S.c])
    print("BC residual", S.bc_residual, "cond", S.bc_cond)
    print("contour eps,kappa,delta (grid) =", eps_n, kap_n, del_n, "Rmax/Rmin/Zmax", ext)
    print("axis x_a =", S.xa, "pb_a =", S.pb_a, "pbxx, pbyy =", S.pbxx, S.pbyy)
    print("psi0 =", S.psi0, "psi_axis =", S.psi_axis, "p_axis =", S.p_axis, "cf =", S.cf)
    print("q_axis (formula) =", S.q_axis(), " q(psiN=1e-6) ray =", S.q_of_psin(1e-6),
          " q_edge =", float(S.q_fit(1.0)), " q_edge ray =", S.q_of_psin(1.0),
          " q(0.5) fit/ray =", float(S.q_fit(0.5)), S.q_of_psin(0.5))
    print("q(0.5) nt=128/512:", S.q_of_psin(0.5, 128), S.q_of_psin(0.5, 512))
    print("Phi_edge (2 pi int q dpsi) =", S.phi_edge, " direct int F/R dA =",
          S.phi_direct(), " ", S.phi_direct(512, 96))
    pv = S.volume_average_p()
    print("<p>_V =", pv, " beta_tor = 2 mu0 <p>/B0^2 =", S.beta_tor(),
          " beta_axis =", 2 * MU0 * S.p_axis / (S.f0 / S.r0) ** 2)
    print("I_tor [A] =", S.toroidal_current())
    import json
    info = dict(c=[float(v) for v in S.c], psi0=S.psi0, psi_axis=S.psi_axis,
                psi_edge=0.0, phi_edge=S.phi_edge, p_axis=S.p_axis, cf=S.cf,
                r_axis=S.r0 * S.xa, z_axis=0.0, q_axis=S.q_axis(),
                q_edge=float(S.q_fit(1.0)), beta_tor_volavg=S.beta_tor(),
                beta_axis=2 * MU0 * S.p_axis / (S.f0 / S.r0) ** 2,
                i_tor_A=S.toroidal_current(), pv_avg=pv, f0=S.f0, r0=S.r0,
                eps=S.eps, kappa=S.kappa, delta=S.delta, A=S.A,
                bc_residual=float(S.bc_residual))
    with open(out("case_info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    print("max |psiN(sample) - psiN(target)| =",
          float(np.max(np.abs(psin - np.repeat(psin_s, 24)))))


if __name__ == "__main__":
    main()
