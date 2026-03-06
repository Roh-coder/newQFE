#!/usr/bin/env python3

import matplotlib.pyplot as plt


def wrapped_points(nx, ny, steps, dx=1, dy=1, x0=0, y0=0):
    pts = []
    for r in range(steps + 1):
        x = (x0 + dx * r) % nx
        y = (y0 + dy * r) % ny
        pts.append((x, y))
    return pts


def unwrapped_points(steps, dx=1, dy=1, x0=0, y0=0):
    return [(x0 + dx * r, y0 + dy * r) for r in range(steps + 1)]


def main():
    nx, ny = 48, 64
    wraps = 2
    steps = wraps * nx

    p_un = unwrapped_points(steps, dx=1, dy=1)
    p_wr = wrapped_points(nx, ny, steps, dx=1, dy=1)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=150)

    # Left: unwrapped trajectory to show linear path in r.
    ax = axes[0]
    xu = [p[0] for p in p_un]
    yu = [p[1] for p in p_un]
    ax.plot(xu, yu, color="#1f77b4", lw=2.0)
    ax.scatter([xu[0], xu[-1]], [yu[0], yu[-1]], c=["green", "red"], s=40, zorder=5)
    ax.text(xu[0], yu[0], " start", color="green", va="bottom", ha="left", fontsize=9)
    ax.text(xu[-1], yu[-1], " end", color="red", va="bottom", ha="left", fontsize=9)
    ax.set_xlabel("unwrapped x = r")
    ax.set_ylabel("unwrapped y = r")
    ax.set_title("x+y trajectory vs separation r (unwrapped)")
    ax.grid(alpha=0.25)

    # Right: wrapped trajectory inside one periodic cell.
    ax = axes[1]
    xw = [p[0] for p in p_wr]
    yw = [p[1] for p in p_wr]

    # Draw periodic box.
    ax.plot([0, nx, nx, 0, 0], [0, 0, ny, ny, 0], color="black", lw=1.2)

    # Draw segments; detect boundary jumps and break line there.
    seg_x = [xw[0]]
    seg_y = [yw[0]]
    for i in range(1, len(xw)):
        dx_jump = abs(xw[i] - xw[i - 1])
        dy_jump = abs(yw[i] - yw[i - 1])
        # Large jump means a periodic wrap happened.
        if dx_jump > 1 or dy_jump > 1:
            ax.plot(seg_x, seg_y, color="#d62728", lw=1.7)
            seg_x = [xw[i]]
            seg_y = [yw[i]]
        else:
            seg_x.append(xw[i])
            seg_y.append(yw[i])
    ax.plot(seg_x, seg_y, color="#d62728", lw=1.7)

    # Mark sampled points every 8 steps for readability.
    mark = list(range(0, len(xw), 8))
    ax.scatter([xw[i] for i in mark], [yw[i] for i in mark], s=12, color="#d62728", alpha=0.9)

    ax.scatter([xw[0], xw[-1]], [yw[0], yw[-1]], c=["green", "red"], s=40, zorder=5)
    ax.text(xw[0], yw[0], " start", color="green", va="bottom", ha="left", fontsize=9)
    ax.text(xw[-1], yw[-1], " end", color="red", va="bottom", ha="left", fontsize=9)

    ax.set_xlim(-1, nx + 1)
    ax.set_ylim(-1, ny + 1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mod 48)")
    ax.set_ylabel("y (mod 64)")
    ax.set_title("x+y trajectory on periodic 48x64 cell (2 wraps in x)")
    ax.grid(alpha=0.2)

    fig.suptitle("Path used by oblique x+y correlation function", fontsize=14)
    fig.tight_layout()

    out_png = "K_from_BC/xplusy_trajectory_L48x64_wrap2.png"
    out_pdf = "K_from_BC/xplusy_trajectory_L48x64_wrap2.pdf"
    fig.savefig(out_png)
    fig.savefig(out_pdf)
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()
