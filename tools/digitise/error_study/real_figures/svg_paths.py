"""Read the stroked paths of a pdftocairo SVG page.

Usage: python3 -I svg_paths.py page.svg   (prints a summary per panel transform)
Imported: paths(svgfile) -> list of dicts with keys transform, stroke, dashed, width, pts (local coordinates).
Only M/L segments are kept (curves and frames in these MATLAB figures are polylines; glyphs are filled, not stroked).
"""
import re, sys

def paths(fn):
    s = open(fn).read()
    out = []
    for attrs in re.findall(r'<path([^>]*)/>', s):
        st = re.search(r'stroke="([^"]*)"', attrs)
        if not st or 'fill="none"' not in attrs:
            continue
        d = re.search(r' d="([^"]*)"', attrs).group(1)
        if re.search(r'[CQAcqa]', d):
            continue
        tr = re.search(r'transform="([^"]*)"', attrs)
        w = re.search(r'stroke-width="([^"]*)"', attrs)
        subs = re.split(r'(?=M )', d.strip())
        for sub in subs:
            pts = [(float(a), float(b)) for a, b in re.findall(r'[ML] ([-\d.e]+) ([-\d.e]+)', sub)]
            if len(pts) < 2:
                continue
            out.append(dict(transform=tr.group(1) if tr else '', stroke=st.group(1),
                            dashed='dasharray' in attrs, width=float(w.group(1)) if w else 1.0, pts=pts))
    return out

if __name__ == '__main__':
    ps = paths(sys.argv[1])
    trs = sorted(set(p['transform'] for p in ps))
    for t in trs:
        print('==', t)
        for p in ps:
            if p['transform'] != t:
                continue
            xs = [a for a, _ in p['pts']]; ys = [b for _, b in p['pts']]
            if len(p['pts']) <= 3 and p['stroke'].startswith('rgb(14.89'):
                continue
            print('  %-42s dash=%d w=%.2f n=%5d x %.2f..%.2f y %.2f..%.2f' % (p['stroke'], p['dashed'], p['width'], len(p['pts']), min(xs), max(xs), min(ys), max(ys)))
