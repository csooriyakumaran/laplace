"""Regenerates docs/fv-cell.png - the node-centred control volume diagram
referenced from README.md §5. Run from the docs/ directory:

    python fv_cell.py
"""
import itertools
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(7.5, 6.5))
ax = fig.add_subplot(111, projection='3d')

corners = list(itertools.product([-1, 1], [-1, 1], [-1, 1]))
for a in corners:
    for b in corners:
        if a < b and sum(1 for i in range(3) if a[i] != b[i]) == 1:
            xs, ys, zs = zip(a, b)
            ax.plot(xs, ys, zs, color='0.3', lw=1.3, zorder=1)

ax.scatter([0], [0], [0], color='black', s=60, zorder=5)
ax.text(-0.95, -0.95, -0.95, '(i, j, k)', fontsize=12, ha='center', zorder=6,
         bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))

faces = [
    ((1, 0, 0),  r'$F_{i+1/2}$'),
    ((-1, 0, 0), r'$F_{i-1/2}$'),
    ((0, 1, 0),  r'$G_{j+1/2}$'),
    ((0, -1, 0), r'$G_{j-1/2}$'),
    ((0, 0, 1),  r'$H_{k+1/2}$'),
    ((0, 0, -1), r'$H_{k-1/2}$'),
]
for (dx, dy, dz), label in faces:
    ax.quiver(0, 0, 0, dx*1.35, dy*1.35, dz*1.35, color='crimson', arrow_length_ratio=0.14, lw=2.2, zorder=4)
    ax.text(dx*1.55, dy*1.55, dz*1.55, label, color='crimson', fontsize=11, ha='center', va='center', zorder=6,
             bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.85))

neighbours = [
    ((2.9, 0, 0), '(i+1, j, k)'), ((-2.9, 0, 0), '(i-1, j, k)'),
    ((0, 2.9, 0), '(i, j+1, k)'), ((0, -2.9, 0), '(i, j-1, k)'),
    ((0, 0, 2.9), '(i, j, k+1)'), ((0, 0, -2.9), '(i, j, k-1)'),
]
for (x, y, z), label in neighbours:
    ax.scatter([x], [y], [z], color='0.55', s=35, zorder=3)
    ax.text(x, y, z+0.22, label, fontsize=9, color='0.35', ha='center', zorder=3)

ax.set_xlabel(r'$\xi\ (i)$', labelpad=-5)
ax.set_ylabel(r'$\eta\ (j)$', labelpad=-5)
ax.set_zlabel(r'$\zeta\ (k)$', labelpad=-5)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])

ax.set_box_aspect((1, 1, 1))
lim = 3.4
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
ax.view_init(elev=16, azim=-50)
ax.set_title('Node-centred control volume (README §5)', pad=0)

plt.tight_layout()
plt.savefig('fv-cell.png', dpi=200, facecolor='white')
