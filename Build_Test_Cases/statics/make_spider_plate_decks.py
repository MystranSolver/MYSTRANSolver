"""Rigid-element timing decks: a 240 x 240 CQUAD4 plate with 900 RBE2 or 900 RBE3 spiders.

The plate is 10 x 10, 241 x 241 grids (grid id = row*241 + col + 1), thickness 0.1, clamped on the edge x = 0. Each spider is
a center grid 0.05 above a 3 x 3 patch of plate grids, with an RBE2 (center dependent, 123456) or an RBE3 (center as reference,
123456, weight 1.0 on 123 of the 9 patch grids), and a unit z force on the center. The patches start at column/row 3 and repeat
every 8 grids: 30 x 30 = 900 spiders, 353,886 G-set DOF and 5,400 M-set DOF. Displacements, SPC forces and MPC forces are printed.

usage: python make_spider_plate_decks.py [OUT_DIR]   (writes spider_plate_rbe2.bdf and spider_plate_rbe3.bdf)
"""
import os
import sys

N = 241                      # grids per side
H = 10.0 / 240               # grid spacing
STEP = 8                     # spider spacing in grids


def plate():
    b = ['SOL 101', 'CEND', 'DISPLACEMENT = ALL', 'SPCFORCE = ALL', 'MPCFORCE = ALL', 'SUBCASE 1', '  LOAD = 1', '  SPC = 1',
         'BEGIN BULK']
    for r in range(N):
        for c in range(N):
            b.append(f'GRID    {r * N + c + 1:<8d}        {c * H:<8.4f}{r * H:<8.4f}0.0')
    for r in range(N - 1):
        for c in range(N - 1):
            g = r * N + c + 1
            b.append(f'CQUAD4  {r * (N - 1) + c + 1:<8d}1       {g:<8d}{g + 1:<8d}{g + N + 1:<8d}{g + N:<8d}')
    b += ['PSHELL  1       1       0.1     1               1', 'MAT1    1       1.0E7           0.3     0.1']
    b += [f'SPC1    1       123456  {r * N + 1:<8d}' for r in range(N)]
    return b


def spiders(kind):
    add, eid = [], 900000
    starts = [s for s in range(3, N - 2, STEP) if s + 2 <= N - 1]
    for r in starts:
        for c in starts:
            g = [(r + i) * N + (c + j) + 1 for i in range(3) for j in range(3)]
            add.append(f'GRID    {eid:<8d}        {(c + 1) * H:<8.4f}{(r + 1) * H:<8.4f}0.05')
            add.append(f'FORCE   1       {eid:<8d}0       1.      0.0     0.0     1.0')
            tag = f'+S{eid - 899999}'
            if kind == 'rbe2':
                add.append(f'RBE2    {eid:<8d}{eid:<8d}123456  ' + ''.join(f'{v:<8d}' for v in g[:5]) + f'{tag:<8s}')
                add.append(f'{tag:<8s}' + ''.join(f'{v:<8d}' for v in g[5:]))
            else:
                add.append(f'RBE3    {eid:<8d}        {eid:<8d}123456  1.      123     ' + ''.join(f'{v:<8d}' for v in g[:2])
                           + f'{tag:<8s}')
                add.append(f'{tag:<8s}' + ''.join(f'{v:<8d}' for v in g[2:]))
            eid += 1
    return add


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    for kind in ('rbe2', 'rbe3'):
        p = os.path.join(out, f'spider_plate_{kind}.bdf')
        open(p, 'w', newline='\n').write('\n'.join(plate() + spiders(kind) + ['PARAM   AUTOSPC YES', 'ENDDATA']) + '\n')
        print('wrote', p)


if __name__ == '__main__':
    main()
