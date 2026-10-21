# Compare every numeric output line of two MYSTRAN F06 files.
# Lines are matched in order after dropping page headers and run-specific lines
# (dates, times, CPU, memory, file names). Each number is compared against a local
# scale: the largest magnitude within 40 lines either side, so round-off noise on
# near-zero terms in a table of large values does not count.
# usage: python -I cmpf06.py a.F06 b.F06 [rtol]
import re, sys

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')
SKIP = re.compile(r'MYSTRAN|PAGE|\d{1,2}/\s?\d{1,2}/\d{4}|\d{1,2}:\s?\d{1,2}:\s?\d{1,2}|CPU|SECONDS|MB |MEGABYTES|ALLOCAT|FILE|ECHO|BEGN|TIME|Total|L\+U', re.I)

def lines(path):
    out = []
    for ln in open(path, errors='replace'):
        if SKIP.search(ln):
            continue
        toks = NUM.findall(ln)
        if len(toks) < 2:
            continue
        try:
            vals = [float(t.replace('D', 'E').replace('d', 'e')) for t in toks]
        except ValueError:
            continue
        out.append((ln.rstrip('\n'), vals))
    return out

def main():
    a, b = lines(sys.argv[1]), lines(sys.argv[2])
    rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5
    if len(a) != len(b):
        print(f'LINE COUNT DIFFERS: {len(a)} vs {len(b)}')
    n = min(len(a), len(b))
    mags = [max((abs(v) for v in a[i][1]), default=0.0) for i in range(n)]
    worst, bad = 0.0, []
    for i in range(n):
        va, vb = a[i][1], b[i][1]
        if len(va) != len(vb):
            bad.append((float('inf'), i)); continue
        scale = max(mags[max(0, i - 40):i + 41]) or 1.0
        r = max(abs(x - y) for x, y in zip(va, vb)) / scale
        worst = max(worst, r)
        if r > rtol:
            bad.append((r, i))
    print(f'numeric lines {n}; worst scaled diff {worst:.2e}; lines over {rtol:g}: {len(bad)}')
    for r, i in sorted(bad, reverse=True)[:5]:
        print(f'  {r:.2e}\n    A: {a[i][0].strip()[:150]}\n    B: {b[i][0].strip()[:150]}')
    return 1 if (bad or len(a) != len(b)) else 0

if __name__ == '__main__':
    sys.exit(main())
