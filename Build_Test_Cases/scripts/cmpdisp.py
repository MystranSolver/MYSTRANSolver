# Compare the grid displacements of two MYSTRAN F06 files (all subcases' tables concatenated by grid).
# usage: python cmpdisp.py a.F06 b.F06
import sys


def disp(path):
    out, on = {}, False
    for ln in open(path, errors='replace'):
        if 'D I S P L A C E M E N T S' in ln:
            on = True
            continue
        if on and ('S P C' in ln or 'F O R C E' in ln or 'S T R E S S' in ln or 'A P P L I E D' in ln):
            on = False
        if on:
            t = ln.split()
            if len(t) == 8 and t[0].isdigit():
                try:
                    out[int(t[0])] = [float(x) for x in t[2:]]
                except ValueError:
                    pass
    return out


a, b = disp(sys.argv[1]), disp(sys.argv[2])
keys = sorted(set(a) & set(b))
mx = max(max(abs(v) for v in a[k]) for k in keys)
d = max(max(abs(x - y) for x, y in zip(a[k], b[k])) for k in keys)
print(f'grids compared {len(keys)} (a {len(a)}, b {len(b)}); max|u| {mx:.6e}; max abs diff {d:.3e}; rel {d / mx:.2e}')
