"""Render the colour fields for the contour-map figures as bare PNGs (axes are drawn in manim)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

D = "data/"
OUT = "images/"
plt.rcParams["font.family"] = "serif"


def field(x, y, Z, name, levels, labels=True, fmt="%.1f", cmap="viridis", vmin=None, vmax=None, logy=False,
          fine=300, label_levels=None):
    """Z indexed [ix, iy]. Saves a borderless image covering exactly x and y ranges."""
    yy = np.log10(y) if logy else y
    interp = RegularGridInterpolator((x, yy), Z, method="linear")
    xf = np.linspace(x.min(), x.max(), fine)
    yf = np.linspace(yy.min(), yy.max(), fine)
    X, Y = np.meshgrid(xf, yf, indexing="ij")
    Zf = interp(np.stack([X, Y], -1))
    fig = plt.figure(figsize=(6, 5), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.contourf(X, Y, Zf, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax, extend="both")
    cs = ax.contour(X, Y, Zf, levels=levels, colors="k", linewidths=0.8)
    if labels:
        ax.clabel(cs, levels=label_levels if label_levels is not None else cs.levels, fmt=fmt, fontsize=13, inline=True)
    ax.set_xlim(xf.min(), xf.max())
    ax.set_ylim(yf.min(), yf.max())
    ax.axis("off")
    fig.savefig(OUT + name + ".png", dpi=200)
    plt.close(fig)


# Fig 1
f1 = np.load(D + "fig1.npz")
field(f1["Ts"], f1["b0"], f1["alpha"].T, "fig1_alpha", np.arange(0.04, 0.19, 0.02), fmt="%.2f", logy=True)
field(f1["Ts"], f1["b0"], np.log10(f1["delta"]).T, "fig1_delta", np.arange(-2.5, 0.51, 0.5), fmt="%.1f", logy=True)

# Fig 6 (shared colour scale 0..1.2)
f6 = np.load(D + "fig6.npz")
lev = np.arange(0.1, 1.25, 0.1)
for k in ["full", "eq26", "eq27"]:
    field(f6["beta"], f6["rs2"], f6[k], "fig6_" + k, lev, vmin=0.1, vmax=1.15)

# Fig 7 (n = 2, 3, 4)
f7 = np.load(D + "fig7.npz")
lev7 = np.arange(0.4, 1.0, 0.05)
for i, n in enumerate([2, 3, 4]):
    field(f7["beta"], f7["rs2"], f7["z830"][:, :, n - 1], f"fig7_n{n}", lev7, fmt="%.2f", vmin=0.1, vmax=1.15)

# Fig 8
f8 = np.load(D + "fig8.npz")
for k in ["thicken", "thin"]:
    field(f8["beta"], f8["rs2"], f8[k], "fig8_" + k, lev, vmin=0.1, vmax=1.15)
print("images done")
