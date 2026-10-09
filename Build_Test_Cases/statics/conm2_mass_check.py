"""Check the CONM2 mass matrix and the grid point weight generator against the exact rigid-mass expressions.

conm2_offset_mass.bdf   CONM2 10. with offset (-3, 0, -2), basic coordinates
conm2_general_mass.bdf  CONM2 10. with offset (1.5, -2, 0.7) and products of inertia I21 = 0.8, I31 = -0.6, I32 = 1.1 in a
                        rotated CID 5, on grid 1 with rotated CP and CD 7; PARAM GRDPNT 0
conm2_cidm1_mass.bdf    the same mass entered with CID = -1 (c.g. coordinates and inertia in basic)

Exact: the inertia tensor J (off-diagonal terms -I21, -I31, -I32 as on the CONM2 entry) and the offset d are rotated to basic,
the mass at the grid is M = m*T'*T + [0 0; 0 J], T = [I, -skew(d)] (the c.g. moves with v + w x d), then rotated to the grid's
CD system. The GPWG inertia about the basic origin is J + m*(|c|^2*I - c*c'), c the c.g. in basic. Both must be met to
rounding; the mass matrix must also be positive semi-definite.

usage: python conm2_mass_check.py --mystran EXE
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

import struct

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def read_op4(path):
    """Binary OUTPUT4, dense columns, as MYSTRAN writes it: {name: matrix}."""
    b, i, mats = open(path, 'rb').read(), 0, {}

    def rec():
        nonlocal i
        n = struct.unpack('<i', b[i:i + 4])[0]
        r = b[i + 4:i + 4 + n]
        i += n + 8
        return r

    while i < len(b):
        h = rec()
        ncol, nrow, form, prec = struct.unpack('<4i', h[:16])
        name = h[16:24].decode().strip()
        a = np.zeros((nrow, ncol))
        while True:
            r = rec()
            j, rbeg, nw = struct.unpack('<3i', r[:12])
            if j > ncol:
                break
            vals = np.frombuffer(r[12:12 + 4 * nw], dtype='<f8')
            a[rbeg - 1:rbeg - 1 + len(vals), j - 1] = vals
        mats[name] = a
    return mats


def cord2r(a, b, c):
    a, b, c = map(np.array, (a, b, c))
    z = (b - a) / np.linalg.norm(b - a)
    x = c - a
    x = x - x.dot(z) * z
    x /= np.linalg.norm(x)
    return a, np.column_stack([x, np.cross(z, x), z])


def skew(d):
    return np.array([[0, -d[2], d[1]], [d[2], 0, -d[0]], [-d[1], d[0], 0]])


def grid_mass(m, d_b, J_b, T_cd):
    T = np.hstack([np.eye(3), -skew(d_b)])
    M = m * T.T @ T
    M[3:, 3:] += J_b
    R = np.zeros((6, 6))
    R[:3, :3] = T_cd
    R[3:, 3:] = T_cd
    return R.T @ M @ R


def run(exe, deck, work):
    os.makedirs(work, exist_ok=True)
    shutil.copy(os.path.join(HERE, deck), work)
    subprocess.run([exe, deck], cwd=work, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
    base = os.path.splitext(deck)[0]
    f06 = open(os.path.join(work, base + '.F06'), errors='replace').read()
    return read_op4(os.path.join(work, base + '.OU1'))['MGG'][:6, :6], f06


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mystran', required=True)
    ap.add_argument('--work', default=os.path.join(HERE, 'conm2_mass_check_out'))
    a = ap.parse_args()
    exe = os.path.abspath(a.mystran)
    ok = True

    def report(name, err, limit):
        nonlocal ok
        good = err <= limit
        ok &= good
        print(f'  {"ok  " if good else "FAIL"} {name}: {err:.2e} (limit {limit:g})')

    G, _ = run(exe, 'conm2_offset_mass.bdf', a.work)
    E = grid_mass(10., np.array([-3., 0., -2.]), np.zeros((3, 3)), np.eye(3))
    print('conm2_offset_mass')
    report('MGG at grid 1 vs exact, of the largest term', abs(G - E).max() / abs(E).max(), 1e-12)
    report('smallest MGG eigenvalue (>= 0), of the largest', max(0.0, -np.linalg.eigvalsh(G).min()) / abs(E).max(), 1e-12)

    G, f06 = run(exe, 'conm2_general_mass.bdf', a.work)
    o5, T5 = cord2r([0, 0, 0], [.3, .4, .866], [1., .2, .1])
    o7, T7 = cord2r([1, 2, 3], [1.2, 1.9, 4.], [2., 2.5, 3.1])
    J5 = np.array([[5., -.8, .6], [-.8, 7., -1.1], [.6, -1.1, 9.]])
    d_b = T5 @ np.array([1.5, -2., .7])
    J_b = T5 @ J5 @ T5.T
    E = grid_mass(10., d_b, J_b, T7)
    print('conm2_general_mass')
    report('MGG at grid 1 vs exact, of the largest term', abs(G - E).max() / abs(E).max(), 1e-12)
    report('smallest MGG eigenvalue (>= 0), of the largest', max(0.0, -np.linalg.eigvalsh(G).min()) / abs(E).max(), 1e-12)
    c = o7 + T7 @ np.array([2., -1., .5]) + d_b
    JO = J_b + 10. * (c.dot(c) * np.eye(3) - np.outer(c, c))
    t = f06.splitlines()
    i = next(k for k, l in enumerate(t) if '6x6 Rigid body mass matrix' in l)
    rows = [r for r in ([float(x) for x in re.findall(r'[-+]?\d\.\d+E[-+]\d+', l)] for l in t[i + 2:i + 9]) if len(r) == 6]
    MO = np.array(rows)
    report('GPWG inertia about the origin vs exact, of the largest term (7-digit print)', abs(MO[3:, 3:] - JO).max() / abs(JO).max(),
           1e-6)

    G1, _ = run(exe, 'conm2_cidm1_mass.bdf', a.work)
    print('conm2_cidm1_mass')
    report('MGG at grid 1 vs exact (CID = -1: c.g. and inertia in basic, 10-digit input), of the largest term',
           abs(G1 - E).max() / abs(E).max(), 1e-8)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
