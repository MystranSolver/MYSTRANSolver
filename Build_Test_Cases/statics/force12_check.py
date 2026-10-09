"""FORCE1, FORCE2, MOMENT1, MOMENT2 against FORCE/MOMENT entries with the same load, and entry names compared exactly.

force12_twins.bdf (written here): a two-bar frame, the loaded grid with a rotated output system (CD 5), four subcases:
  1  FORCE1 (direction G1 to G2) and MOMENT1         2  the same loads as FORCE and MOMENT in basic
  3  FORCE2 (direction (G2 - G1) x (G4 - G3)) and MOMENT2   4  the same loads as FORCE and MOMENT in basic
The directions of subcases 2 and 4 are computed here from the grid coordinates (public NASTRAN definitions). Subcases 1 and 2,
and 3 and 4, must give the same displacements.
prefix_names.bdf (written here): entries whose names only start with an entry MYSTRAN reads (PCOMPG, SUPORT1) must stop the run
with ERROR 1705; earlier versions read them with the shorter entry's layout.
k_shells.bdf (written here): MYSTRAN's own CQUAD4K and CTRIA3K, which the CQUAD4/CTRIA3 readers handle, must still be read;
with --stock, the displacements must equal those of a release build.

usage: python force12_check.py --mystran EXE [--work DIR] [--stock RELEASE_EXE]
"""
import argparse
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
XYZ = {1: (0., 0., 0.), 2: (1., 0., 0.), 3: (2., .5, .3), 10: (.2, .1, .4), 11: (1.3, .9, -.2), 12: (.5, -.4, .7),
       13: (-.3, .6, .1)}


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def d(a, b):
    return np.array(XYZ[b]) - np.array(XYZ[a])


def twins():
    n1 = unit(d(10, 11))                                     # FORCE1 3 100. 10 11
    m1 = unit(d(12, 13))                                     # MOMENT1 3 20. 12 13
    n2 = unit(np.cross(d(10, 11), d(12, 13)))                # FORCE2 3 50. 10 11 12 13
    m2 = unit(np.cross(d(13, 12), d(11, 10)))                # MOMENT2 2 7. 13 12 11 10
    f = lambda v: ','.join(f'{x:.9E}' for x in v)           # noqa: E731
    bulk = ['GRID,1,,0.,0.,0.,,123456', 'GRID,2,,1.,0.,0.', 'GRID,3,,2.,.5,.3,5',
            'CORD2R,5,,0.,0.,0.,.3,.4,.866,+C5', '+C5,1.,.2,.1']
    bulk += [f'GRID,{g},,{x},{y},{z},,123456' for g, (x, y, z) in XYZ.items() if g >= 10]
    bulk += ['CBAR,1,1,1,2,0.,0.,1.', 'CBAR,2,1,2,3,0.,0.,1.', 'PBAR,1,1,1.E-3,2.E-7,3.E-7,4.E-7', 'MAT1,1,7.E10,,.3',
             'FORCE1,1,3,100.,10,11', 'MOMENT1,1,3,20.,12,13',
             f'FORCE,2,3,0,100.,{f(n1)}', f'MOMENT,2,3,0,20.,{f(m1)}',
             'FORCE2,3,3,50.,10,11,12,13', 'MOMENT2,3,2,7.,13,12,11,10',
             f'FORCE,4,3,0,50.,{f(n2)}', f'MOMENT,4,2,0,7.,{f(m2)}']
    case = []
    for k in range(1, 5):
        case += [f'SUBCASE {k}', f'  LOAD = {k}']
    return '\n'.join(['SOL 101', 'CEND', 'TITLE = FORCE1/2 AND MOMENT1/2 TWINS', 'ECHO = NONE', 'DISPLACEMENT(PRINT) = ALL',
                      *case, 'BEGIN BULK', *bulk, 'ENDDATA']) + '\n'


def prefix_names():
    return '\n'.join(['SOL 101', 'CEND', 'TITLE = ENTRY NAMES THAT START WITH A READ ENTRY', 'ECHO = NONE', 'LOAD = 1',
                      'BEGIN BULK', 'GRID,1,,0.,0.,0.,,123456', 'GRID,2,,1.,0.,0.', 'CBAR,1,1,1,2,0.,0.,1.',
                      'PBAR,1,1,1.E-3,2.E-7,3.E-7,4.E-7', 'MAT1,1,7.E10,,.3', 'FORCE,1,2,0,1.,0.,1.,0.',
                      'PCOMPG,7,,,,,,,,+P', '+P,1,1,.1,0.,YES', 'SUPORT1,1,2,3', 'ENDDATA']) + '\n'


def k_shells():
    """CQUAD4K and CTRIA3K (MYSTRAN's own entries, handled by the CQUAD4 and CTRIA3 readers) in a small cantilever plate."""
    b = []
    for j in range(3):
        for i in range(3):
            b.append(f'GRID,{10 * j + i + 1},,{0.5 * i},{0.5 * j},0.' + (',,123456' if i == 0 else ''))
    b += ['CQUAD4K,1,1,1,2,12,11', 'CQUAD4,2,1,11,12,22,21', 'CTRIA3K,3,1,2,3,13', 'CTRIA3,4,1,2,13,12',
          'CQUAD4K,5,1,12,13,23,22', 'PSHELL,1,1,.01,1,,1', 'MAT1,1,7.E10,,.3', 'FORCE,1,23,0,10.,0.,0.,1.',
          'FORCE,1,3,0,5.,0.,0.,1.']
    return '\n'.join(['SOL 101', 'CEND', 'TITLE = CQUAD4K AND CTRIA3K', 'ECHO = NONE', 'LOAD = 1', 'DISPLACEMENT(PRINT) = ALL',
                      'BEGIN BULK', *b, 'ENDDATA']) + '\n'


def run(exe, work, name, text):
    os.makedirs(work, exist_ok=True)
    open(os.path.join(work, name + '.bdf'), 'w', newline='\n').write(text)
    subprocess.run([exe, name + '.bdf'], cwd=work, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
    return open(os.path.join(work, name + '.F06'), errors='replace').read()


def displacements(f06):
    """{subcase: {grid: 6 values}} from the displacement tables."""
    out, sc, on = {}, None, False
    for ln in f06.splitlines():
        m = re.search(r'OUTPUT FOR SUBCASE\s+(\d+)', ln)
        if m:
            sc = int(m.group(1))
        if 'D I S P L A C E M E N T S' in ln:
            on = True
            continue
        if on:
            t = ln.split()
            if len(t) == 8 and t[0].isdigit():
                out.setdefault(sc, {})[int(t[0])] = [float(x) for x in t[2:]]
            elif t and t[0] not in ('GRID', 'COORD') and not t[0].isdigit() and out.get(sc):
                on = False
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mystran', required=True)
    ap.add_argument('--work', default=os.path.join(HERE, 'force12_check_out'))
    ap.add_argument('--stock', default=None, help='a release build to compare the CQUAD4K/CTRIA3K plate with')
    a = ap.parse_args()
    exe = os.path.abspath(a.mystran)
    ok = True
    f06 = run(exe, a.work, 'force12_twins', twins())
    errs = sorted(set(re.findall(r'\*ERROR +\d+', f06)))
    u = displacements(f06)
    for s1, s2, what in ((1, 2, 'FORCE1 + MOMENT1'), (3, 4, 'FORCE2 + MOMENT2')):
        if s1 not in u or s2 not in u:
            print(f'  FAIL {what}: no displacements for subcase {s1} or {s2}; errors {errs}')
            ok = False
            continue
        x = np.array([u[s1][g] for g in sorted(u[s1])])
        y = np.array([u[s2][g] for g in sorted(u[s2])])
        top = abs(y).max()
        err = abs(x - y).max() / top
        good = top > 0 and err <= 1e-5                        # 6-digit F06 print
        ok &= good
        print(f'  {"ok  " if good else "FAIL"} {what} vs FORCE/MOMENT twin: {err:.2e} of the largest displacement ({top:.3e})')
    f06 = run(exe, a.work, 'prefix_names', prefix_names())
    found = sorted(set(re.findall(r'ERROR  1705: MYSTRAN DOES NOT READ BULK DATA ENTRY (\S+)', f06)))
    good = found == ['PCOMPG', 'SUPORT1']
    ok &= good
    print(f'  {"ok  " if good else "FAIL"} ERROR 1705 for {found} (expected PCOMPG, SUPORT1)')
    f06 = run(exe, a.work, 'k_shells', k_shells())
    errs = sorted(set(re.findall(r'\*ERROR +\d+', f06)))
    u = displacements(f06).get(1, {})
    good = not errs and len(u) == 9
    ok &= good
    print(f'  {"ok  " if good else "FAIL"} CQUAD4K and CTRIA3K read and solved: {len(u)} grids, errors {errs}')
    if a.stock:                                           # the same plate with a release build: identical displacements
        ref = displacements(run(os.path.abspath(a.stock), os.path.join(a.work, 'stock'), 'k_shells', k_shells())).get(1, {})
        good = bool(ref) and all(u.get(g) == ref[g] for g in ref)
        ok &= good
        print(f'  {"ok  " if good else "FAIL"} CQUAD4K and CTRIA3K: displacements equal to {os.path.basename(a.stock)}')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
