"""Figure 1: from proximity ligation to the interaction graph.

(a) A promoter and three enhancers assemble in one nucleus; proximity ligation joins two fragments per molecule, and
the capture step keeps only molecules with the baited promoter, so the assembly is read out as three pairs with
exactly two interacting sides each. (b) Transitivity: two measured contacts through a shared fragment imply
proximity of the two unmeasured ends. (c) The pairwise interactions drawn as arcs along the genome become a graph
whose paths (indirect contacts) and densely connected groups (communities) describe proximity beyond the pairs.

Usage: python figure_concept.py <output.pdf>
"""

import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
INK, MUTED, GRID = '#1f1f1e', '#6b6a64', '#d9d8d2'
DASH = (0, (3, 2))


def panel_label(ax, label, title):
    ax.annotate(label, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(0, 4), textcoords='offset points',
                fontweight='bold', fontsize=10, color=INK, va='bottom')
    ax.annotate(title, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(13, 4), textcoords='offset points',
                fontsize=8, color=INK, va='bottom')


def promoter(ax, x, y, name, size=70, label=True):
    ax.scatter(x, y, marker='s', s=size, color=INK, zorder=4)
    if label:
        ax.text(x, y + 0.16, name, ha='center', va='bottom', fontsize=6.5, color=INK, zorder=5)


def other_end(ax, x, y, name, size=60, color='white', label=True):
    ax.scatter(x, y, marker='o', s=size, facecolor=color, edgecolor=INK, lw=1.0, zorder=4)
    if label:
        ax.text(x, y + 0.16, name, ha='center', va='bottom', fontsize=6.5, color=INK, zorder=5)


def molecule(ax, x, y, a, b, kept):
    """A ligation product of two fragments: a and b are 'P' or 'E'."""
    color = INK if kept else GRID
    for dx, kind in ((0.0, a), (0.42, b)):
        if kind == 'P':
            ax.add_patch(mpatches.Rectangle((x + dx, y - 0.09), 0.36, 0.18, color=color, lw=0))
        else:
            ax.add_patch(mpatches.Rectangle((x + dx, y - 0.09), 0.36, 0.18, facecolor='white', edgecolor=color,
                                            lw=1.0))
    ax.plot([x + 0.36, x + 0.42], [y, y], color=color, lw=1.0)


def panel_ligation(ax):
    # nucleus with one promoter-enhancer assembly
    nucleus = mpatches.Ellipse((1.05, 1.5), 2.0, 2.3, facecolor='#f4f3ef', edgecolor=MUTED, lw=0.8)
    ax.add_patch(nucleus)
    t = np.linspace(0, 1, 200)
    fibre_x = 0.25 + 1.6 * t + 0.25 * np.sin(9 * t)
    fibre_y = 0.65 + 1.6 * t ** 1.2 + 0.35 * np.sin(13 * t + 1)
    fibre, = ax.plot(fibre_x, fibre_y, color=MUTED, lw=1.0, zorder=1)
    fibre.set_clip_path(nucleus)
    centre = np.array([1.05, 1.5])
    points = {'P': centre + (0.0, 0.05), 'E1': centre + (-0.32, 0.25), 'E2': centre + (0.3, 0.28),
              'E3': centre + (0.05, -0.3)}
    ax.add_patch(mpatches.Circle(centre + (0.0, 0.07), 0.5, facecolor=AQUA, alpha=0.15, lw=0, zorder=1))
    for (a, b) in (('P', 'E1'), ('P', 'E2'), ('P', 'E3')):
        ax.plot(*zip(points[a], points[b]), color=MUTED, lw=0.6, zorder=2)
    promoter(ax, *points['P'], 'P', label=False)
    for name in ('E1', 'E2', 'E3'):
        other_end(ax, *points[name], name, label=False)
    ax.text(points['P'][0] - 0.22, points['P'][1] - 0.2, 'P', fontsize=6.5, color=INK)
    for name, (dx, dy) in (('E1', (-0.28, 0.0)), ('E2', (0.12, 0.0)), ('E3', (0.12, -0.05))):
        ax.text(points[name][0] + dx, points[name][1] + dy, name, fontsize=6.5, color=INK)
    ax.text(1.05, 0.25, 'one nucleus: promoter\nwith three enhancers', ha='center', va='top', fontsize=6.5,
            color=INK)
    # ligation: every molecule holds two fragments
    ax.annotate('', xy=(2.55, 1.5), xytext=(2.1, 1.5), arrowprops=dict(arrowstyle='->', color=INK, lw=1.0))
    ax.text(2.32, 1.62, 'ligation', ha='center', fontsize=6.5, color=INK)
    pairs = [('P', 'E1', True), ('P', 'E2', True), ('P', 'E3', True), ('E', 'E', False), ('E', 'E', False),
             ('E', 'E', False)]
    names = ['P-E1', 'P-E2', 'P-E3', 'E1-E2', 'E1-E3', 'E2-E3']
    for i, ((a, b, kept), name) in enumerate(zip(pairs, names)):
        y = 2.45 - 0.33 * i
        molecule(ax, 2.7, y, a, b, kept)
        ax.text(3.6, y, name, va='center', fontsize=6.5, color=INK if kept else MUTED)
        ax.text(4.25, y, 'captured' if kept else 'not captured', va='center', fontsize=6.5,
                color=BLUE if kept else MUTED)
    ax.text(2.7, 0.3, 'two fragments per molecule;\nonly promoter molecules captured', fontsize=6.5, color=INK,
            va='bottom')
    ax.set_xlim(0, 5.25)
    ax.set_ylim(-0.3, 2.75)
    ax.set_aspect('equal')
    ax.axis('off')
    panel_label(ax, 'a', 'Proximity ligation records pairs')


def panel_transitivity(ax):
    p, e1, e2 = np.array([1.0, 2.15]), np.array([0.3, 1.0]), np.array([1.7, 1.0])
    ax.plot(*zip(p, e1), color=BLUE, lw=1.8, zorder=2)
    ax.plot(*zip(p, e2), color=BLUE, lw=1.8, zorder=2)
    ax.plot(*zip(e1, e2), color=ORANGE, lw=1.6, ls=DASH, zorder=2)
    promoter(ax, *p, 'v', size=80)
    other_end(ax, *e1, 'u', size=70, label=False)
    other_end(ax, *e2, 'w', size=70, label=False)
    ax.text(e1[0], e1[1] - 0.3, 'u', ha='center', fontsize=7, color=INK)
    ax.text(e2[0], e2[1] - 0.3, 'w', ha='center', fontsize=7, color=INK)
    ax.text(0.38, 1.72, 'measured', rotation=58, ha='center', va='center', fontsize=6.5, color=BLUE)
    ax.text(1.62, 1.72, 'measured', rotation=-58, ha='center', va='center', fontsize=6.5, color=BLUE)
    ax.text(1.0, 0.82, 'inferred', ha='center', fontsize=6.5, color=ORANGE)
    ax.text(1.0, 0.2, 'u ~ v and v ~ w\n$\\Rightarrow$ u near w', ha='center', fontsize=7, color=INK)
    ax.set_xlim(-0.1, 2.1)
    ax.set_ylim(-0.3, 2.75)
    ax.set_aspect('equal')
    ax.axis('off')
    panel_label(ax, 'b', 'Transitivity')


def panel_graph(ax):
    # top: arcs along the genome
    xs = {'E1': 0.2, 'P1': 0.8, 'E2': 1.5, 'E3': 2.0, 'P2': 2.6, 'E4': 3.3, 'E5': 3.9, 'P3': 4.4}
    edges = [('P1', 'E1'), ('P1', 'E2'), ('P1', 'E3'), ('P2', 'E2'), ('P2', 'E3'), ('P2', 'E4'), ('P3', 'E4'),
             ('P3', 'E5'), ('P2', 'P1')]
    base = 2.0

    def arc(x1, x2, height, **kw):
        t = np.linspace(0, np.pi, 60)
        c, h = (x1 + x2) / 2, abs(x2 - x1) / 2
        ax.plot(c - h * np.cos(t), base + height * np.sin(t), **kw)

    ax.plot([0, 4.6], [base, base], color=GRID, lw=0.8, zorder=0)
    for u, v in edges:
        arc(xs[u], xs[v], 0.18 + 0.16 * abs(xs[u] - xs[v]), color=BLUE, lw=1.1)
    for name, x in xs.items():
        if name.startswith('P'):
            promoter(ax, x, base, name, size=40, label=False)
        else:
            other_end(ax, x, base, name, size=34, label=False)
        ax.text(x, base - 0.2, name, ha='center', va='top', fontsize=6, color=INK)
    ax.annotate('', xy=(2.3, 1.38), xytext=(2.3, 1.65), arrowprops=dict(arrowstyle='->', color=INK, lw=1.0))
    # bottom: the same interactions as a graph with communities and indirect contacts
    pos = {'E1': (0.55, 0.55), 'P1': (1.05, 0.85), 'E2': (1.65, 1.05), 'E3': (1.65, 0.5), 'P2': (2.3, 0.8),
           'E4': (3.0, 0.75), 'P3': (3.6, 0.95), 'E5': (3.75, 0.4)}
    ax.add_patch(mpatches.FancyBboxPatch((0.3, 0.25), 2.25, 1.0, boxstyle='round,pad=0.02,rounding_size=0.25',
                                         facecolor=AQUA, alpha=0.13, lw=0))
    ax.add_patch(mpatches.FancyBboxPatch((2.8, 0.18), 1.2, 1.0, boxstyle='round,pad=0.02,rounding_size=0.25',
                                         facecolor=YELLOW, alpha=0.15, lw=0))
    for u, v in (('E1', 'E2'), ('E1', 'E3'), ('E4', 'E2'), ('E5', 'E4')):
        ax.plot(*zip(pos[u], pos[v]), color=ORANGE, lw=1.0, ls=DASH, zorder=2)
    for u, v in edges:
        ax.plot(*zip(pos[u], pos[v]), color=BLUE, lw=1.3, zorder=3)
    for name, (x, y) in pos.items():
        if name.startswith('P'):
            promoter(ax, x, y, name, size=55, label=False)
        else:
            other_end(ax, x, y, name, size=45, label=False)
        dx, dy = {'E1': (-0.12, 0.12), 'P1': (-0.05, 0.13), 'E2': (0.0, 0.12), 'E3': (0.0, -0.22),
                  'P2': (0.05, 0.1), 'E4': (-0.08, 0.13), 'P3': (0.0, 0.12), 'E5': (0.2, -0.08)}[name]
        ax.text(x + dx, y + dy, name, fontsize=6, color=INK, zorder=5, ha='center')
    ax.text(0.35, 0.07, 'community 1', fontsize=6.5, color=AQUA)
    ax.text(2.85, 0.0, 'community 2', fontsize=6.5, color='#b07800')
    ax.set_xlim(-0.05, 4.65)
    ax.set_ylim(-0.3, 2.75)
    ax.set_aspect('equal')
    ax.axis('off')
    panel_label(ax, 'c', 'From arcs to a graph')


def main():
    output = sys.argv[1]
    plt.rcParams.update({'font.size': 8})
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.75), gridspec_kw={'width_ratios': [5.25, 2.2, 4.7]})
    panel_ligation(axes[0])
    panel_transitivity(axes[1])
    panel_graph(axes[2])
    handles = [plt.Line2D([], [], marker='s', color=INK, ls='', label='baited promoter'),
               plt.Line2D([], [], marker='o', mfc='white', mec=INK, ls='', label='other end'),
               plt.Line2D([], [], color=BLUE, lw=1.4, label='measured interaction (edge)'),
               plt.Line2D([], [], color=ORANGE, lw=1.2, ls=DASH, label='indirect contact (two steps)')]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, 0.07))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.08, wspace=0.12)
    fig.savefig(output, bbox_inches='tight', pad_inches=0.02)
    print(output)


if __name__ == '__main__':
    main()
