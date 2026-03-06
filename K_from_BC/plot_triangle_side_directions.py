#!/usr/bin/env python3
"""Plot connected correlators for triangle-side-direction twisted-BC test.

Expected input format (typed .dat):
  type r value err
with connected channels:
  type 4 -> x direction
  type 5 -> oblique side direction (real-y)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def load_connected(path: Path):
    xs, ys = [], []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            t_s, r_s, c_s, e_s = line.split()
            t = int(t_s)
            r = int(r_s)
            c = float(c_s)
            e = float(e_s)
            if t == 4:
                xs.append((r, c, e))
            elif t == 5:
                ys.append((r, c, e))
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--typed', required=True, help='Typed .dat from triangle-side test')
    ap.add_argument('--output', required=True, help='Output image path (png/pdf/svg)')
    ap.add_argument('--title', default='Twisted-BC triangle-side connected correlators')
    args = ap.parse_args()

    typed = Path(args.typed)
    out = Path(args.output)

    xpts, ypts = load_connected(typed)
    if not xpts or not ypts:
        raise SystemExit('Missing connected type-4/type-5 channels in input file.')

    n = min(len(xpts), len(ypts))
    r = [xpts[i][0] for i in range(n)]
    cx = [xpts[i][1] for i in range(n)]
    ex = [xpts[i][2] for i in range(n)]
    cy = [ypts[i][1] for i in range(n)]
    ey = [ypts[i][2] for i in range(n)]

    # Exclude r=0 from difference metric because connected normalization there
    # can dominate and obscure directional matching.
    deltas = []
    for i in range(1, n):
        denom = (ex[i] ** 2 + ey[i] ** 2) ** 0.5
        if denom > 0.0:
            deltas.append((cy[i] - cx[i]) / denom)
    rms_z = (sum(z * z for z in deltas) / len(deltas)) ** 0.5 if deltas else 0.0

    plt.figure(figsize=(8.8, 5.6))
    plt.errorbar(r, cx, yerr=ex, fmt='o-', lw=1.3, ms=3.2, capsize=2.5,
                 label='Connected x (type 4)')
    plt.errorbar(r, cy, yerr=ey, fmt='s-', lw=1.3, ms=3.2, capsize=2.5,
                 label='Connected oblique side / real-y (type 5)')
    plt.axhline(0.0, color='k', lw=0.8, alpha=0.4)
    plt.xlabel('Separation r')
    plt.ylabel('Connected correlator')
    plt.title(f'{args.title}\nRMS z-diff (r>=1): {rms_z:.3f}')
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()

    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=180)
    print(out)


if __name__ == '__main__':
    main()
