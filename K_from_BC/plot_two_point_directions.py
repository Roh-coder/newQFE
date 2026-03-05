#!/usr/bin/env python3
import argparse
from pathlib import Path


def load(path):
    data = {0: [], 1: [], 3: []}
    with open(path) as f:
        for line in f:
            t, r, c, e = line.split()
            t = int(t)
            if t in data:
                data[t].append((int(r), float(c), float(e)))
    return data


def svg_plot(data, title, out_path):
    W, H = 900, 560
    ml, mr, mt, mb = 80, 30, 55, 65
    pw, ph = W - ml - mr, H - mt - mb

    series = []
    labels = {0: 'x direction', 1: 'y direction', 3: 'oblique (x-y) direction'}
    colors = {0: '#1f77b4', 1: '#ff7f0e', 3: '#2ca02c'}

    allx, ally = [], []
    for t in [0, 1, 3]:
        pts = data.get(t, [])
        if not pts:
            continue
        xs = [r for r, _, _ in pts]
        ys = [c for _, c, _ in pts]
        series.append((t, xs, ys))
        allx.extend(xs)
        ally.extend(ys)

    xmin, xmax = min(allx), max(allx)
    ymin, ymax = min(ally), max(ally)
    pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
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

    # axes
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="black"/>')
    parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="black"/>')

    # grid + ticks
    for i in range(6):
        yy = ymin + i * (ymax - ymin) / 5
        ypix = Y(yy)
        parts.append(f'<line x1="{ml}" y1="{ypix:.2f}" x2="{ml+pw}" y2="{ypix:.2f}" stroke="#dddddd"/>')
        parts.append(f'<text x="{ml-10}" y="{ypix+4:.2f}" text-anchor="end" font-size="12" font-family="Arial">{yy:.3f}</text>')

    for i in range(6):
        xx = xmin + i * (xmax - xmin) / 5
        xpix = X(xx)
        parts.append(f'<line x1="{xpix:.2f}" y1="{mt}" x2="{xpix:.2f}" y2="{mt+ph}" stroke="#f0f0f0"/>')
        parts.append(f'<text x="{xpix:.2f}" y="{mt+ph+20}" text-anchor="middle" font-size="12" font-family="Arial">{xx:.1f}</text>')

    parts.append(f'<text x="{ml+pw/2}" y="{H-20}" text-anchor="middle" font-size="14" font-family="Arial">lattice separation bin</text>')
    parts.append(f'<text x="20" y="{mt+ph/2}" transform="rotate(-90 20,{mt+ph/2})" text-anchor="middle" font-size="14" font-family="Arial">C(r)</text>')

    # draw series
    for t, xs, ys in series:
        pts = ' '.join(f'{X(x):.2f},{Y(y):.2f}' for x, y in zip(xs, ys))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{colors[t]}" stroke-width="2"/>')
        for x, y in zip(xs, ys):
            parts.append(f'<circle cx="{X(x):.2f}" cy="{Y(y):.2f}" r="2.5" fill="{colors[t]}"/>')

    # legend
    lx, ly = ml + pw - 210, mt + 15
    parts.append(f'<rect x="{lx}" y="{ly}" width="190" height="75" fill="white" stroke="#cccccc"/>')
    for i, t in enumerate([0, 1, 3]):
        y = ly + 20 + 22 * i
        parts.append(f'<line x1="{lx+10}" y1="{y}" x2="{lx+40}" y2="{y}" stroke="{colors[t]}" stroke-width="3"/>')
        parts.append(f'<text x="{lx+48}" y="{y+4}" font-size="12" font-family="Arial">{labels[t]}</text>')

    parts.append('</svg>')
    Path(out_path).write_text('\n'.join(parts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input')
    ap.add_argument('--output', required=True)
    ap.add_argument('--title', default='Two-point functions by lattice direction')
    args = ap.parse_args()

    data = load(args.input)
    svg_plot(data, args.title, args.output)
    print(args.output)


if __name__ == '__main__':
    main()