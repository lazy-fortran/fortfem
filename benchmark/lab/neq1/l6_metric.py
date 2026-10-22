#!/usr/bin/env python3
"""L6 matched-accuracy vector (FortFEM benchmark/lab/PROBLEMS.md, L6).

Usage: l6_metric.py ref.txt code.txt [--psi-span X]

Both files: columns "R Z psi B_R B_Z B_phi absB p iota" (SI), same rows.  Equal
weights w_i = 1/N, ||x||_w = sqrt(sum_i w_i |x_i|^2) (|.| = Euclidean norm for vector B):
  e_B    = ||B-Bref|| / ||Bref||
  e_p    = ||p-pref|| / max(||pref||, Bref_rms^2/(2 mu0)),  Bref_rms = ||Bref||_w
  e_iota = ||iota-iotaref|| / max(||iotaref||, 1)
  e_flux = ||psi-psiref|| / |psi_edge-psi_axis|
psi_edge-psi_axis is read from the "# psi_axis = ..  psi_edge = .." header of ref.txt
(or --psi-span).  Prints one JSON object with the four numbers and their max.
"""
import argparse
import json
import re

import numpy as np

MU0 = 4e-7 * np.pi
COLUMNS = ("R", "Z", "psi", "B_R", "B_Z", "B_phi", "absB", "p", "iota")


def load(path):
    arr = np.atleast_2d(np.loadtxt(path, comments="#"))
    assert arr.shape[1] == len(COLUMNS), f"{path}: expected {len(COLUMNS)} columns"
    return {k: arr[:, i] for i, k in enumerate(COLUMNS)}


def psi_span(path):
    text = open(path).read()
    m = re.search(r"psi_axis\s*=\s*(\S+)\s+psi_edge\s*=\s*(\S+)", text)
    if not m:
        raise SystemExit(f"{path}: no '# psi_axis = .. psi_edge = ..' header; use --psi-span")
    return abs(float(m.group(2)) - float(m.group(1)))


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


def l6(ref, code, span):
    assert np.allclose(ref["R"], code["R"], rtol=0, atol=1e-12) and np.allclose(
        ref["Z"], code["Z"], rtol=0, atol=1e-12), "sample positions differ"
    n = len(ref["R"])
    bref = np.column_stack([ref["B_R"], ref["B_Z"], ref["B_phi"]])
    b = np.column_stack([code["B_R"], code["B_Z"], code["B_phi"]])
    wn = lambda v: float(np.sqrt(np.mean(np.sum(v**2, axis=1))))  # noqa: E731
    bref_rms = wn(bref)
    out = {
        "e_B": wn(b - bref) / bref_rms,
        "e_p": rms(code["p"] - ref["p"]) / max(rms(ref["p"]), bref_rms**2 / (2 * MU0)),
        "e_iota": rms(code["iota"] - ref["iota"]) / max(rms(ref["iota"]), 1.0),
        "e_flux": rms(code["psi"] - ref["psi"]) / span,
    }
    out["max"] = max(out.values())
    out["n_samples"] = n
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ref")
    ap.add_argument("code")
    ap.add_argument("--psi-span", type=float, default=None)
    a = ap.parse_args()
    span = a.psi_span if a.psi_span else psi_span(a.ref)
    print(json.dumps(l6(load(a.ref), load(a.code), span)))


if __name__ == "__main__":
    main()
