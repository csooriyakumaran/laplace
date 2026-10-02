import sys
from pathlib import Path

import numpy as np
import pandas as pd

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


def add_face_lines(ax, X, Y, Z, **kwargs):
    """Draw every row and every column of a 2D (n, m) node slice as a straight segment"""
    n, m = X.shape
    segments = [np.column_stack([X[a, :], Y[a, :], Z[a, :]]) for a in range(n)]
    segments += [np.column_stack([X[:, b], Y[:, b], Z[:, b]]) for b in range(m)]

    ax.add_collection3d(Line3DCollection(segments, **kwargs))


def main() -> int:

    xs, hs, bs = load_contraction(COORDS)
    X, Y, Z = load_grid(GRID)
    fig, ax = create_fig()

    y_top, z_top = wall_surface(xs, hs, bs)
    x_top = np.broadcast_to(xs[:, None], y_top.shape)
    ax.plot_surface(x_top, y_top, z_top, color='lightsteelblue', alpha=0.4, linewidth=0.4)

    z_side, y_side = wall_surface(xs, bs, hs)
    x_side = np.broadcast_to(xs[:, None], z_side.shape)
    ax.plot_surface(x_side, y_side, z_side, color='lightsteelblue', alpha=0.4, linewidth=0.4)

    cap_kw = dict(color='lightgray', alpha=1.0, linewidth=0)
    ax.plot_surface(X[:, 0, :], Y[:, 0, :], Z[:, 0, :], **cap_kw)  # y-sim plane (j = 0)
    ax.plot_surface(X[:, :, 0], Y[:, :, 0], Z[:, :, 0], **cap_kw)  # z-sim plane (k = 0)
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
