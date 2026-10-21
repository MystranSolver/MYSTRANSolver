"""RBE3 rank decks: the same element in many orientations must give the same verdict.

rbe3_rank_singular.bdf  24 copies of one RBE3 (REFC 123, two weighting grids with 123,
                        reference grid 0.01 off the line through them), each rotated to a
                        different orientation. Rotation about that line is not determined by
                        the weighting grids, so every copy must stop with ERROR 1957.
rbe3_rank_ok.bdf        the same 24 orientations of two well-posed RBE3s: the reference grid
                        on the line (REFC 123), and three non-collinear grids (REFC 123456).
                        The run must finish with no RBE3 error.

Run:  python make_rbe3_rank_decks.py   (writes both decks here)
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def rot(axis, ang):
    x, y, z = axis
    n = math.sqrt(x * x + y * y + z * z)
    x, y, z = x / n, y / n, z / n
    c, s = math.cos(ang), math.sin(ang)
    return [[c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
            [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
            [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)]]


def apply(r, p):
    return [sum(r[i][k] * p[k] for k in range(3)) for i in range(3)]


def f8(v):
    for p in range(7, 0, -1):
        s = f'{v:.{p}g}'
        if 'e' in s:
            mant, exp = s.split('e')
            s = mant + ('' if '.' in mant else '.') + f'E{int(exp):+d}'
        elif '.' not in s:
            s += '.'
        if len(s) <= 8:
            return s.ljust(8)
    raise ValueError(v)


def grid(gid, p, ps=''):
    return f'GRID    {gid:<8d}        {f8(p[0])}{f8(p[1])}{f8(p[2])}        {ps}'.rstrip()


ORIENTS = [rot(a, ang) for a in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 1, 1), (0.3, -0.8, 0.5))
           for ang in (0.0, 0.4, 1.1, math.pi / 2)]


def deck(title, body):
    head = ['SOL 101', 'CEND', f'TITLE = {title}', 'ECHO = NONE', 'SUBCASE 1', '  LOAD = 1', '  SPC = 1',
            '  DISPLACEMENT = ALL', 'BEGIN BULK']
    tail = ['MAT1    1       7.0E10          0.33    2700.', 'PROD    1       1       1.0E-4',
            'PARAM   AUTOSPC YES', 'ENDDATA', '']
    return '\n'.join(head + body + tail)


def singular():
    b = []
    for k, r in enumerate(ORIENTS):
        base = 1000 * (k + 1)
        pts = {1: [0.0, 0.0, 0.0], 2: [0.1, 0.0, 0.0], 3: [0.05, 0.01, 0.0]}   # 3 = reference, off the line 1-2
        for i, p in pts.items():
            b.append(grid(base + i, apply(r, p), '123456' if i < 3 else ''))
        b.append(f'RBE3    {base + 1:<8d}        {base + 3:<8d}123     1.      123     {base + 1:<8d}{base + 2:<8d}')
        b.append(f'FORCE   1       {base + 3:<8d}0       1.      1.      1.      1.')
    b.append('SPC1    1       123456  999999')
    b.append(grid(999998, [8.0, 9.0, 9.0], '123456'))      # MYSTRAN needs at least one element
    b.append(grid(999999, [9.0, 9.0, 9.0], '123456'))
    b.append('CROD    1       1       999998  999999')
    return deck('RBE3 RANK: TWO GRIDS, REFERENCE OFF THEIR LINE, 24 ORIENTATIONS (ALL MUST FAIL 1957)', b)


def ok():
    b = []
    for k, r in enumerate(ORIENTS):
        base = 1000 * (k + 1)
        pts = {1: [0.0, 0.0, 0.0], 2: [0.1, 0.0, 0.0], 3: [0.04, 0.0, 0.0],    # 3 = reference on the line 1-2
               4: [0.0, 0.1, 0.0], 5: [0.03, 0.03, 0.0]}                      # 1, 2, 4 non-collinear; 5 = reference
        for i, p in pts.items():
            b.append(grid(base + i, apply(r, p), '123456' if i in (1, 2, 4) else ''))
        b.append(f'RBE3    {base + 1:<8d}        {base + 3:<8d}123     1.      123     {base + 1:<8d}{base + 2:<8d}')
        cont = f'+R{k:<6d}'                             # third grid on a continuation (field 10 is the marker)
        b.append(f'RBE3    {base + 2:<8d}        {base + 5:<8d}123456  1.      123     {base + 1:<8d}{base + 2:<8d}{cont}')
        b.append(f'{cont}{base + 4:<8d}')
        b.append(f'FORCE   1       {base + 3:<8d}0       1.      1.      1.      1.')
        b.append(f'FORCE   1       {base + 5:<8d}0       1.      1.      1.      1.')
    b.append('SPC1    1       123456  999999')
    b.append(grid(999998, [8.0, 9.0, 9.0], '123456'))      # MYSTRAN needs at least one element
    b.append(grid(999999, [9.0, 9.0, 9.0], '123456'))
    b.append('CROD    1       1       999998  999999')
    return deck('RBE3 RANK: WELL-POSED ELEMENTS, 24 ORIENTATIONS (NO RBE3 ERROR)', b)


if __name__ == '__main__':
    for name, make in (('rbe3_rank_singular.bdf', singular), ('rbe3_rank_ok.bdf', ok)):
        with open(os.path.join(HERE, name), 'w', newline='\n') as f:
            f.write(make())
        print('wrote', name)
