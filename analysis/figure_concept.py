"""Figure 1: from proximity ligation to the interaction graph.

(a) A promoter and three enhancers assemble in one nucleus; proximity ligation joins two fragments per molecule, and
the capture step keeps only molecules with the baited promoter, so the assembly is read out as three pairs with
exactly two interacting sides each. (b) Indirect contact: two measured contacts through a shared fragment connect
the two ends, which were never measured as a pair. (c) The pairwise interactions drawn as arcs along the genome become a graph of the measured
interactions; in the same graph with the measured interactions in grey, all pairs at two steps (indirect contacts)
and the communities are shown.

Usage: python figure_concept.py <output prefix>   (writes <prefix>a.pdf, <prefix>b.pdf, <prefix>c.pdf)
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
    """Panel title; the panel letter is set by LaTeX when the panels are assembled."""
    ax.annotate(title, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(0, 4), textcoords='offset points',
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
    # nucleus with one promoter-enhancer assembly and the DNA strand
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
    ax.text(1.05, 0.25, 'promoter with\nthree enhancers', ha='center', va='top', fontsize=6.5, color=INK)
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


def panel_indirect(ax):
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
    ax.text(1.0, 0.82, 'indirect', ha='center', fontsize=6.5, color=ORANGE)
    ax.text(1.0, 0.2, 'u - v and v - w measured\n$\\Rightarrow$ u - w indirect contact', ha='center', fontsize=6.5, color=INK)
    ax.set_xlim(-0.1, 2.1)
    ax.set_ylim(-0.3, 2.75)
    ax.set_aspect('equal')
    ax.axis('off')
    panel_label(ax, 'b', 'Indirect contact')


EDGES = [('P1', 'E1'), ('P1', 'E2'), ('P1', 'E3'), ('P2', 'E2'), ('P2', 'E3'), ('P2', 'E4'), ('P3', 'E4'),
         ('P3', 'E5'), ('P2', 'P1')]
# planar layout: P1 and P2 on a midline with their shared partners E2 above and E3 below, so that no two edges cross
LAYOUT = {'E1': (0.5, 0.7), 'P1': (1.05, 0.7), 'E2': (1.68, 1.12), 'E3': (1.68, 0.28), 'P2': (2.3, 0.7),
          'E4': (3.0, 0.7), 'P3': (3.62, 1.05), 'E5': (3.62, 0.35)}
LABELS = {'E1': (-0.12, -0.25), 'P1': (0.0, -0.25), 'E2': (0.0, 0.11), 'E3': (0.0, -0.24), 'P2': (0.0, -0.25),
          'E4': (0.0, -0.25), 'P3': (0.0, 0.11), 'E5': (0.0, -0.24)}
# bow of each indirect contact (perpendicular offset of the curve's control point; 0 = straight)
BOWS = {('E1', 'E2'): 0.25, ('E1', 'E3'): -0.25, ('E1', 'P2'): 1.3, ('E2', 'E3'): 0.0, ('E4', 'P1'): 1.3,
        ('E2', 'E4'): 0.25, ('E3', 'E4'): -0.25, ('P2', 'P3'): 0.0, ('E4', 'E5'): 0.0}


def two_step_pairs(edges):
    neighbours = {}
    for u, v in edges:
        neighbours.setdefault(u, set()).add(v)
        neighbours.setdefault(v, set()).add(u)
    pairs = set()
    for v, around in neighbours.items():
        for u in around:
            for w in around:
                if u < w and w not in neighbours[u]:
                    pairs.add((u, w))
    return sorted(pairs)


def draw_graph(ax, ox, mode):
    """mode 'measured': edges in blue; 'indirect': edges grey, two-step contacts in orange, communities shaded."""
    pos = {n: (x + ox, y) for n, (x, y) in LAYOUT.items()}
    if mode == 'indirect':
        ax.add_patch(mpatches.FancyBboxPatch((ox + 0.3, -0.02), 2.2, 1.44, boxstyle='round,pad=0.02,rounding_size=0.25',
                                             facecolor=AQUA, alpha=0.13, lw=0))
        ax.add_patch(mpatches.FancyBboxPatch((ox + 2.78, -0.02), 1.12, 1.44,
                                             boxstyle='round,pad=0.02,rounding_size=0.25',
                                             facecolor=YELLOW, alpha=0.15, lw=0))
        ax.text(ox + 0.35, -0.25, 'community 1', fontsize=6.5, color=AQUA)
        ax.text(ox + 2.8, -0.25, 'community 2', fontsize=6.5, color='#b07800')
        t = np.linspace(0, 1, 60)[:, None]
        for u, w in two_step_pairs(EDGES):
            bow = BOWS.get((u, w), BOWS.get((w, u), 0.0))
            (x1, y1), (x2, y2) = pos[u], pos[w]
            length = np.hypot(x2 - x1, y2 - y1)
            control = np.array([(x1 + x2) / 2 - bow * (y2 - y1) / length, (y1 + y2) / 2 + bow * (x2 - x1) / length])
            points = (1 - t) ** 2 * np.array([x1, y1]) + 2 * (1 - t) * t * control + t ** 2 * np.array([x2, y2])
            ax.plot(points[:, 0], points[:, 1], color=ORANGE, lw=1.1, ls=DASH, zorder=3)
    for u, v in EDGES:
        ax.plot(*zip(pos[u], pos[v]), color=BLUE if mode == 'measured' else GRID, lw=1.3, zorder=2)
    for name, (x, y) in pos.items():
        if name.startswith('P'):
            promoter(ax, x, y, name, size=55, label=False)
        else:
            other_end(ax, x, y, name, size=45, label=False)
        dx, dy = LABELS[name]
        ax.text(x + dx, y + dy, name, fontsize=6, color=INK, zorder=5, ha='center')


def panel_graph(ax):
    # left: arcs along the genome
    xs = {'E1': 0.2, 'P1': 0.8, 'E2': 1.5, 'E3': 2.0, 'P2': 2.6, 'E4': 3.3, 'E5': 3.9, 'P3': 4.4}
    base = 0.45

    def arc(x1, x2, height, **kw):
        t = np.linspace(0, np.pi, 60)
        c, h = (x1 + x2) / 2, abs(x2 - x1) / 2
        ax.plot(c - h * np.cos(t), base + height * np.sin(t), **kw)

    ax.plot([0, 4.6], [base, base], color=GRID, lw=0.8, zorder=0)
    for u, v in EDGES:
        arc(xs[u], xs[v], 0.18 + 0.16 * abs(xs[u] - xs[v]), color=BLUE, lw=1.1)
    for name, x in xs.items():
        if name.startswith('P'):
            promoter(ax, x, base, name, size=40, label=False)
        else:
            other_end(ax, x, base, name, size=34, label=False)
        ax.text(x, base - 0.2, name, ha='center', va='top', fontsize=6, color=INK)
    ax.text(2.3, -0.16, 'pairwise interactions along the genome', ha='center', fontsize=6.5, color=MUTED)
    for x in (4.85, 9.35):
        ax.annotate('', xy=(x + 0.45, 0.7), xytext=(x, 0.7), arrowprops=dict(arrowstyle='->', color=INK, lw=1.0))
    draw_graph(ax, 5.2, 'measured')
    ax.text(5.2 + 2.2, -0.16, 'graph of the measured interactions', ha='center', fontsize=6.5, color=MUTED)
    draw_graph(ax, 9.7, 'indirect')
    ax.text(9.7 + 2.2, 1.55, 'indirect contacts and communities', ha='center', fontsize=6.5, color=MUTED)
    ax.set_xlim(-0.05, 13.7)
    ax.set_ylim(-0.35, 1.7)
    ax.set_aspect('equal')
    ax.axis('off')
    panel_label(ax, 'c', 'From arcs to a graph')


def save(fig, path):
    fig.savefig(path, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    print(path)


def main():
    prefix = sys.argv[1]
    plt.rcParams.update({'font.size': 8})
    for letter, panel, size in (('a', panel_ligation, (4.6, 2.75)), ('b', panel_indirect, (1.95, 2.75))):
        fig, ax = plt.subplots(figsize=size)
        panel(ax)
        save(fig, '{}{}.pdf'.format(prefix, letter))
    fig, ax = plt.subplots(figsize=(7.2, 1.75))
    panel_graph(ax)
    handles = [plt.Line2D([], [], marker='s', color=INK, ls='', label='baited promoter'),
               plt.Line2D([], [], marker='o', mfc='white', mec=INK, ls='', label='other end'),
               plt.Line2D([], [], color=BLUE, lw=1.4, label='measured interaction (edge)'),
               plt.Line2D([], [], color=ORANGE, lw=1.2, ls=DASH, label='indirect contact (two steps)')]
    ax.legend(handles=handles, loc='upper center', ncol=4, frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, -0.02))
    save(fig, prefix + 'c.pdf')


if __name__ == '__main__':
    main()
