#!/usr/bin/env python3
import argparse
from pathlib import Path


def load(path):
    data = {0: [], 1: [], 3: [], 4: [], 5: []}
    with open(path) as f:
        for line in f:
            t, r, c, e = line.split()
            t = int(t)
            if t in data:
                data[t].append((int(r), float(c), float(e)))
    return data


def svg_plot(data, title, out_path, length_by_type=None, caption=''):
    W, H = 900, 590
    ml, mr, mt, mb = 80, 30, 55, 95
    pw, ph = W - ml - mr, H - mt - mb

    series = []
    labels = {
        0: 'x direction',
        1: 'y direction',
        3: 'oblique (x-y) direction',
        4: 'connected x (jackknife)',
        5: 'connected y (jackknife)',
    }
    colors = {0: '#1f77b4', 1: '#ff7f0e', 3: '#2ca02c', 4: '#1f77b4', 5: '#ff7f0e'}

    allx, ally = [], []
    if data.get(4) or data.get(5):
        draw_order = [4, 5, 3]
    else:
        draw_order = [0, 1, 3]

    for t in draw_order:
        pts = data.get(t, [])
        if not pts:
            continue
        if length_by_type is not None and t in length_by_type and length_by_type[t] > 0:
            L = float(length_by_type[t])
            xs = [r / L for r, _, _ in pts]
        else:
            xs = [r for r, _, _ in pts]
        ys = [c for _, c, _ in pts]
        es = [e for _, _, e in pts]
        series.append((t, xs, ys, es))
        allx.extend(xs)
        ally.extend(ys)
        for (_x, y, e) in pts:
            ally.append(y - e)
            ally.append(y + e)

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

    x_label = 'lattice separation bin'
    if length_by_type is not None:
        x_label = 'fractional separation r/L'
    parts.append(f'<text x="{ml+pw/2}" y="{H-44}" text-anchor="middle" font-size="14" font-family="Arial">{x_label}</text>')
    if caption:
        parts.append(f'<text x="{ml+pw/2}" y="{H-20}" text-anchor="middle" font-size="12" font-family="Arial">{caption}</text>')
    parts.append(f'<text x="20" y="{mt+ph/2}" transform="rotate(-90 20,{mt+ph/2})" text-anchor="middle" font-size="14" font-family="Arial">C(r)</text>')

    # draw series
    cap = 4.0
    for t, xs, ys, es in series:
        # Draw vertical error bars and caps first so line/markers remain visible.
        for x, y, e in zip(xs, ys, es):
            xpix = X(x)
            ylo = Y(y - e)
            yhi = Y(y + e)
            parts.append(
                f'<line x1="{xpix:.2f}" y1="{ylo:.2f}" x2="{xpix:.2f}" y2="{yhi:.2f}" '
                f'stroke="{colors[t]}" stroke-width="1.2" opacity="0.85"/>'
            )
            parts.append(
                f'<line x1="{xpix-cap:.2f}" y1="{ylo:.2f}" x2="{xpix+cap:.2f}" y2="{ylo:.2f}" '
                f'stroke="{colors[t]}" stroke-width="1.2" opacity="0.85"/>'
            )
            parts.append(
                f'<line x1="{xpix-cap:.2f}" y1="{yhi:.2f}" x2="{xpix+cap:.2f}" y2="{yhi:.2f}" '
                f'stroke="{colors[t]}" stroke-width="1.2" opacity="0.85"/>'
            )

        pts = ' '.join(f'{X(x):.2f},{Y(y):.2f}' for x, y in zip(xs, ys))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{colors[t]}" stroke-width="2"/>')
        for x, y in zip(xs, ys):
            parts.append(f'<circle cx="{X(x):.2f}" cy="{Y(y):.2f}" r="2.5" fill="{colors[t]}"/>')

    # legend
    lx, ly = ml + pw - 210, mt + 15
    legend_types = [t for t in draw_order if data.get(t)]
    legend_h = max(30, 20 + 22 * len(legend_types))
    parts.append(f'<rect x="{lx}" y="{ly}" width="190" height="{legend_h}" fill="white" stroke="#cccccc"/>')
    for i, t in enumerate(legend_types):
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
    ap.add_argument('--normalize_fraction', action='store_true',
                    help='Plot x-axis as fractional separation r/L (requires --nx/--ny for types 0/1).')
    ap.add_argument('--nx', type=int, default=0,
                    help='Lx for type 0 (x-direction) when using --normalize_fraction.')
    ap.add_argument('--ny', type=int, default=0,
                    help='Ly for type 1 (y-direction) when using --normalize_fraction.')
    ap.add_argument('--caption', default='',
                    help='Optional caption line rendered below the x-axis label.')
    args = ap.parse_args()

    data = load(args.input)
    length_by_type = None
    if args.normalize_fraction:
        length_by_type = {}
        if args.nx > 0:
            length_by_type[0] = args.nx
            length_by_type[3] = args.nx
            length_by_type[4] = args.nx
        if args.ny > 0:
            length_by_type[1] = args.ny
            length_by_type[5] = args.ny
    svg_plot(data, args.title, args.output, length_by_type=length_by_type, caption=args.caption)
    print(args.output)


if __name__ == '__main__':
    main()