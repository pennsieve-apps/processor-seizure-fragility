"""256-bit reference norms for the cases written by make_cases.R (python-flint ball arithmetic).

For each <case>_A.bin (int32 n, int32 nSearch, A column-major, Re and Im of the trial
eigenvalues, all little-endian doubles) and each trial eigenvalue lambda,
M = (A - lambda I)^-1 is computed with arb/acb at 256 bits from the exact doubles, and for each
row r of M the minimum-norm perturbation of EZFragility's fragilityRow (Delta with
[Im r; Re r] Delta = (0, -1)) has |Delta|^2 = |Im r|^2 / (|Im r|^2 |Re r|^2 - (Im r . Re r)^2),
evaluated in the same precision, where its cancellation is harmless. Writes <case>_norms.bin
(n x nSearch doubles, column-major) and prints the largest relative radius of the balls, a
rigorous bound on the reference's own error.

usage: python3 tests/reference/reference.py tests/testthat/reference
"""
import os
import struct
import sys

import numpy as np
from flint import acb, acb_mat, arb, ctx

ctx.prec = 256


def norms(A, lam):
    n = A.shape[0]
    out = np.full((n, len(lam)), np.nan)
    rad = 0.0
    for k, (lr, li) in enumerate(lam):
        l = acb(arb(float(lr)), arb(float(li)))
        M = acb_mat(n, n)
        for i in range(n):
            for j in range(n):
                M[i, j] = acb(arb(float(A[i, j])))
            M[i, i] = M[i, i] - l
        Minv = M.inv()
        for i in range(n):
            ii = arb(0); rr = arb(0); ir = arb(0)
            for j in range(n):
                e = Minv[i, j]
                ii += e.imag * e.imag; rr += e.real * e.real; ir += e.imag * e.real
            d = ii * rr - ir * ir
            if d > 0:
                v = (ii / d).sqrt()
                out[i, k] = float(v.mid())
                rad = max(rad, float(v.rad()) / abs(float(v.mid())))
    return out, rad


def main(folder):
    for f in sorted(os.listdir(folder)):
        if not f.endswith("_A.bin"):
            continue
        with open(os.path.join(folder, f), "rb") as fh:
            n, ns = struct.unpack("<ii", fh.read(8))
            A = np.frombuffer(fh.read(8 * n * n), "<f8").reshape((n, n), order="F")
            lv = np.frombuffer(fh.read(16 * ns), "<f8")
        V, rad = norms(A, list(zip(lv[:ns], lv[ns:])))
        V.astype("<f8").ravel(order="F").tofile(os.path.join(folder, f.replace("_A.bin", "_norms.bin")))
        print("%s: %d x %d norms, max relative radius %.1e" % (f, n, ns, rad))


if __name__ == "__main__":
    main(sys.argv[1])
