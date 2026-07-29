"""
Shared plotting style for the matched-null audit.

Usage in any notebook:
    from audit_style import *
    # Then use FIG_SINGLE, FIG_DOUBLE, OUTLIER_COLOR, etc.
    # apply_style() is called automatically on import.

For publication-ready figures, save as PDF:
    fig.savefig('figs/my_figure.pdf', bbox_inches='tight')
"""

import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Repository-relative paths. Scripts live in analysis/; everything they read
# and write is inside this repository.
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "outputs"
CACHE = REPO / "cache"
for _d in (OUT, CACHE):
    _d.mkdir(parents=True, exist_ok=True)
BASE = REPO  # legacy alias


# ── Colour palette ─────────────────────────────────────────────────────────
OUTLIER_COLOR  = '#E05A4F'   # tomato-ish, distinguishable in B&W
CONTROL_COLOR = '#4682B4'   # steelblue
NEUTRAL_COLOR = '#888888'   # grey for background/reference elements
ACCENT_COLOR  = '#2CA02C'   # green for third category if needed

PALETTE = {'outlier': OUTLIER_COLOR, 'control': CONTROL_COLOR}
PALETTE_NEURONAL = {
    'outlier_neuronal':  OUTLIER_COLOR,
    'control_neuronal': CONTROL_COLOR,
    'non_neuronal':     NEUTRAL_COLOR,
}

# ── Standard figure sizes (inches) ────────────────────────────────────────
# Adjust to match your LaTeX template column widths
FIG_SINGLE = (3.5, 3.0)     # single-column figure
FIG_DOUBLE = (7.0, 4.0)     # double-column figure
FIG_FULL   = (7.0, 7.0)     # full-page panel figure
FIG_WIDE   = (10.0, 4.0)    # wide figure (presentations)
FIG_PANEL  = (14.0, 11.0)   # multi-panel (4+ subplots)

# ── Output directories ────────────────────────────────────────────────────
FIGS_DIR = Path('figs')
DATA_DIR = Path('data')

def apply_style():
    """Apply the project's plotting style globally.

    Uses seaborn 'paper' context with custom overrides tuned for
    Nature Methods single- and double-column figures.
    """
    sns.set(style='white', context='paper', font_scale=1.1)
    plt.rcParams.update({
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'pdf.fonttype': 42,       # TrueType in PDF (editable in Illustrator)
        'ps.fonttype': 42,
        'font.size': 9,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.bottom': True,       # Nature Methods requires tick marks
        'ytick.left': True,
        'legend.frameon': False,
        'lines.linewidth': 0.8,
    })
    FIGS_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)


def save_fig(fig, name, formats=('pdf', 'png')):
    """Save figure in multiple formats.

    Args:
        fig: matplotlib Figure object
        name: filename stem (no extension), e.g. 'gf_01_pca_outliers'
        formats: tuple of format strings to save
    """
    for fmt in formats:
        path = FIGS_DIR / f'{name}.{fmt}'
        fig.savefig(path, dpi=300, bbox_inches='tight')
    fmts = ', '.join(f'figs/{name}.{f}' for f in formats)
    print(f'Saved: {fmts}')


def repel_labels(ax, texts, x_coords, y_coords, max_iter=100):
    """Simple label repulsion to reduce text overlap.

    Works by iteratively nudging overlapping text annotations apart.
    A lightweight alternative to adjustText when that package is unavailable.

    Args:
        ax: matplotlib Axes
        texts: list of matplotlib Text objects (from ax.annotate or ax.text)
        x_coords: array of original x positions (data coords)
        y_coords: array of original y positions (data coords)
        max_iter: number of repulsion iterations
    """
    import numpy as np

    if len(texts) < 2:
        return

    fig = ax.get_figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for _ in range(max_iter):
        moved = False
        bboxes = [t.get_window_extent(renderer=renderer) for t in texts]
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                bi, bj = bboxes[i], bboxes[j]
                # Check overlap
                if bi.overlaps(bj):
                    # Push apart vertically in display coords
                    overlap_y = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
                    shift = (overlap_y / 2 + 2)  # pixels
                    # Shift in data coordinates
                    inv = ax.transData.inverted()
                    _, dy = inv.transform((0, shift)) - inv.transform((0, 0))
                    pos_i = texts[i].get_position()
                    pos_j = texts[j].get_position()
                    texts[i].set_position((pos_i[0], pos_i[1] + dy))
                    texts[j].set_position((pos_j[0], pos_j[1] - dy))
                    moved = True
        if not moved:
            break
        # Refresh bboxes
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()


# Apply style on import
apply_style()
