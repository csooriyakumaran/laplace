import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import aero

DATA_DIR = Path(__file__).parent.parent / 'data'
COORDS = DATA_DIR / 'contraction.csv'
GRID = DATA_DIR / 'output.csv'


def load_contraction(path):
    df = pd.read_csv(path)
    return df['xs'].to_numpy(), df['hs'].to_numpy(), df['bs'].to_numpy()


def load_grid(path):
    df = pd.read_csv(path)
    ni, nj, nk = df['i'].max() + 1, df['j'].max() + 1, df['k'].max() + 1

    X = np.empty((ni, nj, nk))
    Y = np.empty((ni, nj, nk))
    Z = np.empty((ni, nj, nk))

    i, j, k = df['i'].to_numpy(), df['j'].to_numpy(), df['k'].to_numpy()

    X[i, j, k] = df['x'].to_numpy()
    Y[i, j, k] = df['y'].to_numpy()
    Z[i, j, k] = df['z'].to_numpy()

    return X, Y, Z


def load_velocity(path):
    """u, v, w columns written by laplace_output (application/c/main.c), same
    (ni, nj, nk) node layout as load_grid. Returns None if the CSV predates
    that output (grid-only, no velocity columns)."""
    df = pd.read_csv(path)
    if not {'u', 'v', 'w'}.issubset(df.columns):
        return None

    ni, nj, nk = df['i'].max() + 1, df['j'].max() + 1, df['k'].max() + 1
    U = np.empty((ni, nj, nk))
    V = np.empty((ni, nj, nk))
    W = np.empty((ni, nj, nk))

    i, j, k = df['i'].to_numpy(), df['j'].to_numpy(), df['k'].to_numpy()

    U[i, j, k] = df['u'].to_numpy()
    V[i, j, k] = df['v'].to_numpy()
    W[i, j, k] = df['w'].to_numpy()

    return U, V, W


def create_fig() -> tuple[aero.plot.AeroPlot3D, aero.plot.Axes]:
    spec = aero.plot.Spec3D()
    spec.title = 'Laplace 1/4 Domain Grid'
    spec.xlabels = ['$x$']
    spec.ylabels = ['$y$']
    spec.zlabels = ['$z$']
    spec.zlims = []
    spec.view_kw = dict(elev=35, azim=25, roll=0, vertical_axis='z')

    fig: aero.plot.AeroPlot3D = aero.plot.figure(spec)
    ax: aero.plot.Axes = fig.axis(row=0, col=0)

    return fig, ax


def wall_surface(xs, along_vals, extent_vals, n_sweep=60):
    """Ruled surface: 'along' = along_val(x), sweep over [0, extent_vals], in the cross direction"""
    s = np.linspace(0.0, 1.0, n_sweep)
    cross = np.outer(extent_vals, s)
    along = np.broadcast_to(along_vals[:, None], cross.shape)
    return along, cross


def add_velocity_on_plane(fig, ax, X, Y, Z, U, V, W, plane='k0', stride=1, cmap_name='viridis'):
    """Illustrate velocity on one symmetry plane: a speed-colored surface
    plus a direction-only quiver (uniform arrow length, colored by the same
    speed scale would double-encode magnitude, so arrows just show direction).

    plane: 'k0' = z-sym plane (k = 0, spanwise centreline); 'j0' = y-sym
    plane (j = 0, vertical centreline).
    """
    if plane == 'k0':
        Xs, Ys, Zs = X[:, :, 0], Y[:, :, 0], Z[:, :, 0]
        Us, Vs, Ws = U[:, :, 0], V[:, :, 0], W[:, :, 0]
    elif plane == 'j0':
        Xs, Ys, Zs = X[:, 0, :], Y[:, 0, :], Z[:, 0, :]
        Us, Vs, Ws = U[:, 0, :], V[:, 0, :], W[:, 0, :]
    else:
        raise ValueError("plane must be 'k0' or 'j0'")

    speed = np.sqrt(Us**2 + Vs**2 + Ws**2)

    norm = Normalize(vmin=speed.min(), vmax=speed.max())
    cmap = plt.colormaps[cmap_name]
    ax.plot_surface(Xs, Ys, Zs, facecolors=cmap(norm(speed)), shade=False, linewidth=0, antialiased=False)

    sl = (slice(None, None, stride), slice(None, None, stride))
    # ax.quiver(Xs[sl], Ys[sl], Zs[sl], Us[sl], Vs[sl], Ws[sl],
    #           length=0.15 * max(X.max() - X.min(), Y.max() - Y.min()), normalize=True, color='k', linewidth=0.6)

    bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6, pad=0.1)
    bar.ax.set_ylabel('$|V|$')


def add_face_lines(ax, X, Y, Z, **kwargs):
    """Draw every row and every column of a 2D (n, m) node slice as a straight segment"""
    n, m = X.shape
    segments = [np.column_stack([X[a, :], Y[a, :], Z[a, :]]) for a in range(n)]
    segments += [np.column_stack([X[:, b], Y[:, b], Z[:, b]]) for b in range(m)]

    ax.add_collection3d(Line3DCollection(segments, **kwargs))


def plot_axis_velocity(x, y, z, u, v, w):
    mask = np.where((y == 0) & (z == 0))
    fig = aero.plot.figure(aero.plot.Spec2D())

    U = np.sqrt(u * u + v * v + w * w)

    fig.add_scatter(x[:, 0, -1], u[:, 0, -1], label='wall (u)')
    fig.add_scatter(x[:, 0, 0], u[:, 0, 0], label='axis (u)')
    fig.add_scatter(x[:, 0, 0], U[:, 0, 0], label='axis (U)')
    fig.add_legend()


def plot_inlet_velocity(x, y, z, u, v, w):
    fig = aero.plot.figure(aero.plot.Spec2D())
    fig.add_scatter(u[0, :, 0], y[0, :, 0], label='inlet')
    fig.add_legend()


def plot_outlet_velocity(x, y, z, u, v, w):
    fig = aero.plot.figure(aero.plot.Spec2D())
    fig.add_scatter(u[-1, :, 0], y[-1, :, 0], label='outlet')
    fig.add_legend()


def main() -> int:

    xs, hs, bs = load_contraction(COORDS)
    X, Y, Z = load_grid(GRID)
    velocity = load_velocity(GRID)
    fig, ax = create_fig()

    y_top, z_top = wall_surface(xs, hs, bs)
    x_top = np.broadcast_to(xs[:, None], y_top.shape)
    ax.plot_surface(x_top, y_top, z_top, color='lightsteelblue', alpha=0.4, linewidth=0.4)

    z_side, y_side = wall_surface(xs, bs, hs)
    x_side = np.broadcast_to(xs[:, None], z_side.shape)
    ax.plot_surface(x_side, y_side, z_side, color='lightsteelblue', alpha=0.4, linewidth=0.4)

    cap_kw = dict(color='lightgray', alpha=1.0, linewidth=0)
    ax.plot_surface(X[:, 0, :], Y[:, 0, :], Z[:, 0, :], **cap_kw)  # y-sim plane (j = 0)
    if velocity is None:
        ax.plot_surface(X[:, :, 0], Y[:, :, 0], Z[:, :, 0], **cap_kw)  # z-sim plane (k = 0)
    else:
        U, V, W = velocity
        add_velocity_on_plane(fig, ax, X, Y, Z, U, V, W, plane='k0')  # z-sym plane (k = 0), colored by speed
        plot_axis_velocity(X, Y, Z, U, V, W)
        plot_inlet_velocity(X, Y, Z, U, V, W)
        plot_outlet_velocity(X, Y, Z, U, V, W)

    ax.plot_surface(X[0, :, :], Y[0, :, :], Z[0, :, :], **cap_kw)  # inlet
    ax.plot_surface(X[-1, :, :], Y[-1, :, :], Z[-1, :, :], **cap_kw)  # outlet

    line_kw = dict(colors='magenta', linewidths=0.4)
    # add_face_lines(ax, X[:, -1, :], Y[:, -1, :], Z[:, -1, :], **line_kw)  # top wall (j = nj - 1)
    # add_face_lines(ax, X[:, :, -1], Y[:, :, -1], Z[:, :, -1], **line_kw)  # side wall (k = nk - 1)
    add_face_lines(ax, X[:, 0, :], Y[:, 0, :], Z[:, 0, :], **line_kw)  # y-sym plane (j = 0)
    add_face_lines(ax, X[:, :, 0], Y[:, :, 0], Z[:, :, 0], **line_kw)  # z-sym plane (k = 0)
    add_face_lines(ax, X[0, :, :], Y[0, :, :], Z[0, :, :], **line_kw)  # inlet (i = 0)
    add_face_lines(ax, X[-1, :, :], Y[-1, :, :], Z[-1, :, :], **line_kw)  # outlet (i = ni - 1)

    # ax.axis('off')
    fig.set_axes_equal()
    aero.plot.show()

    return 0


if __name__ == '__main__':
    sys.exit(main())
