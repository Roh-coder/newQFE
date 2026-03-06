#!/usr/bin/env python3
import argparse
from pathlib import Path


def load_typed(path):
    data = {0: [], 1: []}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            t, r, c, e = line.split()
            t = int(t)
            if t in data:
                data[t].append((int(r), float(c), float(e)))
    return data


def svg_plot(base, test, title, out_path, nx, ny, normalize_fraction=False):
    W, H = 980, 620
    ml, mr, mt, mb = 90, 40, 60, 75
    pw, ph = W - ml - mr, H - mt - mb

    series = []
    # (label, type, source, color)
    spec = [
        ('baseline x', 0, base, '#1f77b4'),
        ('baseline y', 1, base, '#ff7f0e'),
        ('h/b=3/4 x', 0, test, '#2ca02c'),
        ('h/b=3/4 y', 1, test, '#d62728'),
    ]

    allx, ally = [], []
    for label, t, src, color in spec:
        pts = src.get(t, [])
        if not pts:
            continue
        L = nx if t == 0 else ny
        if normalize_fraction:
            xs = [r / float(L) for r, _, _ in pts]
        else:
            xs = [r for r, _, _ in pts]
        ys = [c for _, c, _ in pts]
        es = [e for _, _, e in pts]
        series.append((label, xs, ys, es, color))
        allx.extend(xs)
        ally.extend(ys)
        for y, e in zip(ys, es):
            ally.append(y - e)
            ally.append(y + e)

    xmin, xmax = min(allx), max(allx)
    ymin, ymax = min(ally), max(ally)
    pad = 0.06 * (ymax - ymin if ymax > ymin else 1.0)
    ymin -= pad
    ymax += pad

    def X(x):
        return ml + (x - xmin) / (xmax - xmin if xmax > xmin else 1.0) * pw

    def Y(y):
        return mt + (ymax - y) / (ymax - ymin if ymax > ymin else 1.0) * ph

    parts = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">')
    parts.append('<rect x="0" y="0" width="100%" height="100%" fill="white"/>')
    parts.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-size="20" font-family="Arial">{title}</text>')

    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="black"/>')
    parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="black"/>')

    for i in range(6):
        yy = ymin + i * (ymax - ymin) / 5
        ypix = Y(yy)
        parts.append(f'<line x1="{ml}" y1="{ypix:.2f}" x2="{ml+pw}" y2="{ypix:.2f}" stroke="#e0e0e0"/>')
        parts.append(f'<text x="{ml-10}" y="{ypix+4:.2f}" text-anchor="end" font-size="12" font-family="Arial">{yy:.3f}</text>')

    for i in range(6):
        xx = xmin + i * (xmax - xmin) / 5
        xpix = X(xx)
        parts.append(f'<line x1="{xpix:.2f}" y1="{mt}" x2="{xpix:.2f}" y2="{mt+ph}" stroke="#f0f0f0"/>')
        lbl = f'{xx:.2f}' if normalize_fraction else f'{xx:.1f}'
        parts.append(f'<text x="{xpix:.2f}" y="{mt+ph+20}" text-anchor="middle" font-size="12" font-family="Arial">{lbl}</text>')

    xlab = 'fractional separation r/L' if normalize_fraction else 'separation r'
    parts.append(f'<text x="{ml+pw/2}" y="{H-20}" text-anchor="middle" font-size="14" font-family="Arial">{xlab}</text>')
    parts.append(f'<text x="20" y="{mt+ph/2}" transform="rotate(-90 20,{mt+ph/2})" text-anchor="middle" font-size="14" font-family="Arial">C_conn(r)</text>')

    cap = 3.5
    for label, xs, ys, es, color in series:
        for x, y, e in zip(xs, ys, es):
            xpix = X(x)
            ylo = Y(y - e)
            yhi = Y(y + e)
            parts.append(f'<line x1="{xpix:.2f}" y1="{ylo:.2f}" x2="{xpix:.2f}" y2="{yhi:.2f}" stroke="{color}" stroke-width="1" opacity="0.85"/>')
            parts.append(f'<line x1="{xpix-cap:.2f}" y1="{ylo:.2f}" x2="{xpix+cap:.2f}" y2="{ylo:.2f}" stroke="{color}" stroke-width="1" opacity="0.85"/>')
            parts.append(f'<line x1="{xpix-cap:.2f}" y1="{yhi:.2f}" x2="{xpix+cap:.2f}" y2="{yhi:.2f}" stroke="{color}" stroke-width="1" opacity="0.85"/>')

        pts = ' '.join(f'{X(x):.2f},{Y(y):.2f}' for x, y in zip(xs, ys))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2"/>')

    lx, ly = ml + pw - 220, mt + 12
    parts.append(f'<rect x="{lx}" y="{ly}" width="205" height="95" fill="white" stroke="#cccccc"/>')
    for i, (label, _xs, _ys, _es, color) in enumerate(series):
        y = ly + 18 + 20 * i
        parts.append(f'<line x1="{lx+10}" y1="{y}" x2="{lx+40}" y2="{y}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{lx+48}" y="{y+4}" font-size="12" font-family="Arial">{label}</text>')

    parts.append('</svg>')
    Path(out_path).write_text('\n'.join(parts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--baseline', required=True)
    ap.add_argument('--test', required=True)
    ap.add_argument('--nx', type=int, required=True)
    ap.add_argument('--ny', type=int, required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--title', default='Connected two-point comparison')
    ap.add_argument('--normalize_fraction', action='store_true')
    args = ap.parse_args()

    base = load_typed(args.baseline)
    test = load_typed(args.test)
    svg_plot(base, test, args.title, args.output, args.nx, args.ny,
             normalize_fraction=args.normalize_fraction)
    print(args.output)


if __name__ == '__main__':
    main()
