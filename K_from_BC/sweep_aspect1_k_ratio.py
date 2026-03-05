#!/usr/bin/env python3
import argparse
import math
import struct
import subprocess
import zlib
from pathlib import Path


THETA_ASPECT1_ISOSCELES = math.atan(2.0)


def load_xy(dat_path):
    x, y = [], []
    with open(dat_path) as f:
        for line in f:
            t, r, c, _e = line.split()
            t = int(t)
            rv = int(r)
            cv = float(c)
            ev = float(_e)
            if t == 0:
                x.append((rv, cv, ev))
            elif t == 1:
                y.append((rv, cv, ev))
    return x, y


def write_xy_dat(path, xpts, ypts):
    n = min(len(xpts), len(ypts))
    with open(path, 'w') as f:
        f.write('# r corr_x err_x corr_y err_y\n')
        for i in range(n):
            rx, cx, ex = xpts[i]
            ry, cy, ey = ypts[i]
            if rx != ry:
                raise ValueError(f'mismatched bins at index {i}: x has {rx}, y has {ry}')
            f.write(f'{rx} {cx:.16e} {ex:.16e} {cy:.16e} {ey:.16e}\n')


def write_png(path, rgb, width, height):
    def chunk(tag, data):
        return (struct.pack('!I', len(data)) + tag + data +
                struct.pack('!I', zlib.crc32(tag + data) & 0xFFFFFFFF))

    raw = bytearray()
    stride = width * 3
    for row in range(height):
        raw.append(0)
        start = row * stride
        raw.extend(rgb[start:start + stride])

    png = bytearray(b'\x89PNG\r\n\x1a\n')
    png += chunk(b'IHDR', struct.pack('!IIBBBBB', width, height, 8, 2, 0, 0, 0))
    png += chunk(b'IDAT', zlib.compress(bytes(raw), level=9))
    png += chunk(b'IEND', b'')
    Path(path).write_bytes(png)


def _write_simple_pdf(path, width, height, content_stream):
    objects = []
    objects.append(b'<< /Type /Catalog /Pages 2 0 R >>')
    objects.append(b'<< /Type /Pages /Kids [3 0 R] /Count 1 >>')
    page = (
        f'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} {height}] '
        f'/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>'
    )
    objects.append(page.encode())
    objects.append(b'<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>')
    stream = content_stream.encode('ascii')
    objects.append(b'<< /Length ' + str(len(stream)).encode() + b' >>\nstream\n' + stream + b'\nendstream')

    out = bytearray(b'%PDF-1.4\n')
    offsets = [0]
    for i, obj in enumerate(objects, start=1):
        offsets.append(len(out))
        out.extend(f'{i} 0 obj\n'.encode())
        out.extend(obj)
        out.extend(b'\nendobj\n')

    xref_start = len(out)
    out.extend(f'xref\n0 {len(objects) + 1}\n'.encode())
    out.extend(b'0000000000 65535 f \n')
    for off in offsets[1:]:
        out.extend(f'{off:010d} 00000 n \n'.encode())

    trailer = f'trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_start}\n%%EOF\n'
    out.extend(trailer.encode())
    Path(path).write_bytes(out)


def render_xy_pdf(xpts, ypts, title, out_path):
    width, height = 1000, 650
    ml, mr, mt, mb = 90, 40, 70, 80
    pw, ph = width - ml - mr, height - mt - mb

    all_x = [r for r, _, _ in xpts] + [r for r, _, _ in ypts]
    all_y = [c for _, c, _ in xpts] + [c for _, c, _ in ypts]
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
    ymin -= pad
    ymax += pad

    def X(x):
        return ml + (x - xmin) / (xmax - xmin if xmax > xmin else 1.0) * pw

    def Y(y):
        return mt + (ymax - y) / (ymax - ymin if ymax > ymin else 1.0) * ph

    def py(y_top):
        return height - y_top

    cmd = []
    cmd.append('1 1 1 rg 0 0 1000 650 re f')

    cmd.append('0.88 0.88 0.88 RG 0.8 w')
    for i in range(6):
        yy = ymin + i * (ymax - ymin) / 5
        ypix = Y(yy)
        cmd.append(f'{ml:.2f} {py(ypix):.2f} m {ml + pw:.2f} {py(ypix):.2f} l S')
    for i in range(6):
        xx = xmin + i * (xmax - xmin) / 5
        xpix = X(xx)
        cmd.append(f'{xpix:.2f} {py(mt):.2f} m {xpix:.2f} {py(mt + ph):.2f} l S')

    cmd.append('0 0 0 RG 1.2 w')
    cmd.append(f'{ml:.2f} {py(mt + ph):.2f} m {ml + pw:.2f} {py(mt + ph):.2f} l S')
    cmd.append(f'{ml:.2f} {py(mt):.2f} m {ml:.2f} {py(mt + ph):.2f} l S')

    def draw_series(pts, rgb):
        r, g, b = rgb
        cmd.append(f'{r:.3f} {g:.3f} {b:.3f} RG 1.6 w')
        pix = [(X(px), Y(pyv)) for px, pyv, _ in pts]
        if len(pix) >= 2:
            x0, y0 = pix[0]
            cmd.append(f'{x0:.2f} {py(y0):.2f} m')
            for x1, y1 in pix[1:]:
                cmd.append(f'{x1:.2f} {py(y1):.2f} l')
            cmd.append('S')

    draw_series(xpts, (31 / 255.0, 119 / 255.0, 180 / 255.0))
    draw_series(ypts, (255 / 255.0, 127 / 255.0, 14 / 255.0))

    safe_title = title.replace('\\', '\\\\').replace('(', '\\(').replace(')', '\\)')
    cmd.append('0 0 0 rg BT /F1 16 Tf 90 620 Td (' + safe_title + ') Tj ET')

    _write_simple_pdf(out_path, width, height, '\n'.join(cmd))


def render_xy_png(xpts, ypts, title, out_path):
    width, height = 1000, 650
    ml, mr, mt, mb = 90, 40, 70, 80
    pw, ph = width - ml - mr, height - mt - mb

    all_x = [r for r, _, _ in xpts] + [r for r, _, _ in ypts]
    all_y = [c for _, c, _ in xpts] + [c for _, c, _ in ypts]
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
    ymin -= pad
    ymax += pad

    img = bytearray([255] * (width * height * 3))

    def set_px(x, y, color):
        if 0 <= x < width and 0 <= y < height:
            i = (y * width + x) * 3
            img[i:i + 3] = bytes(color)

    def line(x0, y0, x1, y1, color):
        dx = abs(x1 - x0)
        sx = 1 if x0 < x1 else -1
        dy = -abs(y1 - y0)
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        while True:
            set_px(x0, y0, color)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    def circle(cx, cy, rad, color):
        for yy in range(cy - rad, cy + rad + 1):
            for xx in range(cx - rad, cx + rad + 1):
                if (xx - cx) ** 2 + (yy - cy) ** 2 <= rad * rad:
                    set_px(xx, yy, color)

    def X(x):
        return int(round(ml + (x - xmin) / (xmax - xmin if xmax > xmin else 1.0) * pw))

    def Y(y):
        return int(round(mt + (ymax - y) / (ymax - ymin if ymax > ymin else 1.0) * ph))

    black = (0, 0, 0)
    grid = (225, 225, 225)
    blue = (31, 119, 180)
    orange = (255, 127, 14)

    line(ml, mt + ph, ml + pw, mt + ph, black)
    line(ml, mt, ml, mt + ph, black)

    for i in range(6):
        yy = ymin + i * (ymax - ymin) / 5
        ypix = Y(yy)
        line(ml, ypix, ml + pw, ypix, grid)

    for i in range(6):
        xx = xmin + i * (xmax - xmin) / 5
        xpix = X(xx)
        line(xpix, mt, xpix, mt + ph, grid)

    for pts, color in ((xpts, blue), (ypts, orange)):
        pixels = [(X(r), Y(c)) for r, c, _ in pts]
        for (x0, y0), (x1, y1) in zip(pixels[:-1], pixels[1:]):
            line(x0, y0, x1, y1, color)
        for x, y in pixels:
            circle(x, y, 3, color)

    write_png(out_path, img, width, height)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', default='/tmp/ising_kbc_aspect1')
    ap.add_argument('--nx', type=int, default=24)
    ap.add_argument('--ny', type=int, default=24)
    ap.add_argument('--k1', type=float, default=1.0)
    ap.add_argument('--ratios', default='0.70,0.85,1.00,1.15,1.30')
    ap.add_argument('--n_therm', type=int, default=300)
    ap.add_argument('--n_traj', type=int, default=2200)
    ap.add_argument('--n_skip', type=int, default=10)
    ap.add_argument('--n_wolff', type=int, default=3)
    ap.add_argument('--n_metropolis', type=int, default=2)
    ap.add_argument('--out_dir', default='K_from_BC/aspect1_scan_pdf')
    ap.add_argument('--xy_dat_dir', default='K_from_BC/aspect1_scan_xy_dat',
                    help='Directory for extracted x/y correlator .dat files')
    ap.add_argument('--formats', default='pdf',
                    help='Comma-separated list: pdf,png')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    xy_dat_dir = Path(args.xy_dat_dir)
    xy_dat_dir.mkdir(parents=True, exist_ok=True)

    ratios = [float(r.strip()) for r in args.ratios.split(',') if r.strip()]
    formats = {f.strip().lower() for f in args.formats.split(',') if f.strip()}

    print('ratio,mid_bin,corr_x_mid,corr_y_mid,y_over_x,xy_dat,outputs')
    for ratio in ratios:
        k2 = args.k1 * ratio
        k3 = args.k1 * ratio
        cmd = [
            args.exe,
            '--n_x', str(args.nx), '--n_y', str(args.ny),
            '--theta1', str(THETA_ASPECT1_ISOSCELES),
            '--theta2', str(THETA_ASPECT1_ISOSCELES),
            '--k1', str(args.k1), '--k2', str(k2), '--k3', str(k3),
            '--n_therm', str(args.n_therm), '--n_traj', str(args.n_traj),
            '--n_skip', str(args.n_skip), '--n_wolff', str(args.n_wolff),
            '--n_metropolis', str(args.n_metropolis),
        ]
        run = subprocess.run(cmd, check=True, text=True, capture_output=True)

        dat = Path(
            f'ising_flat_crit_k_from_bc/{args.nx}_{args.ny}_t1_{THETA_ASPECT1_ISOSCELES:.3f}_t2_{THETA_ASPECT1_ISOSCELES:.3f}_k_{args.k1:.3f}_{k2:.3f}_{k3:.3f}.dat'
        )
        xpts, ypts = load_xy(dat)

        midpoint = len(xpts) // 2
        cx_mid = xpts[midpoint][1]
        cy_mid = ypts[midpoint][1]
        ratio_mid = cy_mid / cx_mid if abs(cx_mid) > 0.0 else float('nan')

        xy_dat = xy_dat_dir / f'aspect1_k1_{args.k1:.2f}_ratio_{ratio:.2f}_xy.dat'
        write_xy_dat(xy_dat, xpts, ypts)

        title = f'Aspect=1 scan: k1={args.k1:.2f}, k2/k1=k3/k1={ratio:.2f}'
        outputs = []
        if 'pdf' in formats:
            pdf = out_dir / f'aspect1_k1_{args.k1:.2f}_ratio_{ratio:.2f}.pdf'
            render_xy_pdf(xpts, ypts, title, pdf)
            outputs.append(str(pdf))
        if 'png' in formats:
            png = out_dir / f'aspect1_k1_{args.k1:.2f}_ratio_{ratio:.2f}.png'
            render_xy_png(xpts, ypts, title, png)
            outputs.append(str(png))

        log = out_dir / f'run_ratio_{ratio:.2f}.log'
        log.write_text(run.stdout)

        print(f'{ratio:.3f},{midpoint},{cx_mid:.6e},{cy_mid:.6e},{ratio_mid:.6f},{xy_dat},{";".join(outputs)}')


if __name__ == '__main__':
    main()