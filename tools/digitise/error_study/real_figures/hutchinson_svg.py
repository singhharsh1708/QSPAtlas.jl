"""Exact curve coordinates of Hutchinson 2016 Figures 2 and 6 from the vector figures in the first
author's DPhil thesis (Oxford ORA uuid:0c1034ee-ff8e-4ed4-b168-f559f93f9142): thesis Figure 3.2
(PDF page 85) is paper Figure 2, thesis Figure 3.14 (PDF page 111) is paper Figure 6.

Prepare the SVGs (not part of the repository) with
  pdftocairo -svg -f 85 -l 85 hutchinson_thesis.pdf p85.svg
  pdftocairo -svg -f 111 -l 111 hutchinson_thesis.pdf p111.svg
Run:
  python3 -I hutchinson_svg.py <folder with p85.svg and p111.svg> <output folder>

Calibration: every tick mark is a 2-point path in the panel's own coordinates. x ticks are the
vertical ticks on the bottom frame, y ticks the horizontal ticks on the left frame. Tick values
are the printed labels (first value and step per panel, below); the map is a least-squares line
through all ticks and the largest tick residual is reported (in data units).
Writes <output>/<atlas id>.csv (x,y in data units, every vertex of the drawn polyline) and
<output>/calibration.txt.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from svg_paths import paths

BLUE = 'rgb(0%, 44.709778%, 74.119568%)'; ORANGE = 'rgb(85.099792%, 32.548523%, 9.799194%)'
YELLOW = 'rgb(92.939758%, 69.40918%, 12.548828%)'; PURPLE = 'rgb(49.409485%, 18.429565%, 55.688477%)'
GREEN = 'rgb(46.669006%, 67.449951%, 18.818665%)'; BLACK = 'rgb(0%, 0%, 0%)'; MAGENTA = 'rgb(100%, 0%, 100%)'
ORDER = [BLUE, ORANGE, YELLOW, PURPLE, GREEN]   # legend order of thesis Fig 3.2 (c)-(f)

# (svg, transform offsets, frame left, frame bottom, x first, x step, y first (bottom tick), y step, curves {id: colour})
F32 = 'p85.svg'; F314 = 'p111.svg'
PANELS = [
    (F32, (116.962745, 85.391115), 83.0, 355.0, 40, 5, 50, 50, {'fig2_a_T': BLUE}),
    (F32, (299.22507, 85.391115), 83.0, 355.0, 40, 5, 100, 100, {'fig2_b_V': BLUE}),
    (F32, (116.962745, 222.323572), 83.0, 355.0, 40, 5, 50, 50, dict(zip(['fig2_c_T_d04', 'fig2_c_T_d06', 'fig2_c_T_d08', 'fig2_c_T_d10', 'fig2_c_T_d12'], ORDER))),
    (F32, (299.22507, 222.323572), 83.0, 355.0, 40, 5, 50, 50, dict(zip(['fig2_d_V_d04', 'fig2_d_V_d06', 'fig2_d_V_d08', 'fig2_d_V_d10', 'fig2_d_V_d12'], ORDER))),
    (F32, (116.962745, 359.257028), 83.0, 355.0, 40, 5, 50, 50, dict(zip(['fig2_e_T_N2', 'fig2_e_T_N3', 'fig2_e_T_N6', 'fig2_e_T_N11', 'fig2_e_T_N21'], ORDER))),
    (F32, (299.22507, 359.257028), 83.0, 355.0, 40, 5, 50, 50, dict(zip(['fig2_f_V_N2', 'fig2_f_V_N3', 'fig2_f_V_N6', 'fig2_f_V_N11', 'fig2_f_V_N21'], ORDER))),
    (F314, (136.782639, 85.391714), 83.0, 275.0, 30, 10, 0, 40, {'fig6_a_T': BLACK, 'fig6_a_V': MAGENTA}),
    (F314, (299.22507, 85.391714), 83.0, 269.0, 35, 5, 0, 40, {'fig6_b_T': BLACK, 'fig6_b_V': MAGENTA}),
    (F314, (136.782639, 190.113347), 83.0, 269.0, 35, 5, 0, 40, {'fig6_c_T': BLACK, 'fig6_c_V': MAGENTA}),
    (F314, (299.22507, 190.113347), 83.0, 269.0, 35, 5, 0, 40, {'fig6_d_T': BLACK, 'fig6_d_V': MAGENTA}),
]

def lsq(p, v):
    n = len(p); mp = sum(p) / n; mv = sum(v) / n
    b = sum((a - mp) * (c - mv) for a, c in zip(p, v)) / sum((a - mp) ** 2 for a in p)
    a = mv - b * mp
    return a, b, max(abs(a + b * x - y) for x, y in zip(p, v))

def main(src, out):
    os.makedirs(out, exist_ok=True)
    cache = {}
    log = open(os.path.join(out, 'calibration.txt'), 'w')
    for svg, off, left, bottom, x0, dx, y0, dy, curves in PANELS:
        ps = cache.setdefault(svg, paths(os.path.join(src, svg)))
        mine = [p for p in ps if p['transform'].endswith('%s, %s)' % off)]
        assert mine, (svg, off)
        two = [p['pts'] for p in mine if len(p['pts']) == 2]
        def tlen(q): return max(abs(q[0][0] - q[1][0]), abs(q[0][1] - q[1][1]))
        xt = sorted({round(q[0][0], 4) for q in two if 3 < tlen(q) < 8 and abs(q[0][0] - q[1][0]) < 1e-6 and min(abs(q[0][1] - bottom), abs(q[1][1] - bottom)) < 0.05})
        yt = sorted({round(q[0][1], 4) for q in two if 3 < tlen(q) < 8 and abs(q[0][1] - q[1][1]) < 1e-6 and min(abs(q[0][0] - left), abs(q[1][0] - left)) < 0.05}, reverse=True)
        ax, bx, rx = lsq(xt, [x0 + i * dx for i in range(len(xt))])
        ay, by, ry = lsq(yt, [y0 + i * dy for i in range(len(yt))])
        log.write('%s %s: %d x ticks %g..%g (max residual %.2e), %d y ticks %g..%g (max residual %.2e)\n' % (svg, off, len(xt), x0, x0 + dx * (len(xt) - 1), rx, len(yt), y0, y0 + dy * (len(yt) - 1), ry))
        for cid, col in curves.items():
            segs = [p['pts'] for p in mine if p['stroke'] == col and not p['dashed'] and len(p['pts']) >= 10]
            assert segs, cid
            pts = sorted({pt for s in segs for pt in s})
            with open(os.path.join(out, cid + '.csv'), 'w') as f:
                f.write('x,y\n')
                for X, Y in pts:
                    f.write('%.6g,%.6g\n' % (ax + bx * X, ay + by * Y))
            log.write('  %s: %d vertices from %d path(s)\n' % (cid, len(pts), len(segs)))
    log.close()

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
