"""Shared publication styling for the complete November 12 figure set."""
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gwas_submission_matplotlib_20261112')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPOSITORY = Path(__file__).resolve().parents[2]
ROOT = REPOSITORY / 'figures/submission_2026-11-12'
DATA_ROOT = REPOSITORY / 'derived_data/submission_2026-11-12'
METHOD_ORDER = ['BOLT-LMM', 'REGENIE_LOCO', 'REGENIE_WO_LOCO', 'fastGWA',
                'FarmCPU_NYC', 'FarmCPU_Kansas', 'BLINK_NYC', 'BLINK_Kansas',
                'SLEMM', 'GCTA_MLMA']
METHOD_LABELS = {
    'BOLT-LMM': 'BOLT-LMM', 'REGENIE_LOCO': 'REGENIE_LOCO',
    'REGENIE_WO_LOCO': 'REGENIE_woLOCO', 'fastGWA': 'fastGWA',
    'FarmCPU_NYC': 'FarmCPU (NYC)', 'FarmCPU_Kansas': 'FarmCPU (Kansas)',
    'BLINK_NYC': 'BLINK (NYC)', 'BLINK_Kansas': 'BLINK (Kansas)',
    'SLEMM': 'SLEMM', 'GCTA_MLMA': 'GCTA-MLMA',
}
METHOD_COLORS = {
    'BOLT-LMM': '#0072B2', 'REGENIE_LOCO': '#D55E00',
    'REGENIE_WO_LOCO': '#E69F00', 'fastGWA': '#009E73',
    'FarmCPU_NYC': '#D81B60', 'FarmCPU_Kansas': '#004D40',
    'BLINK_NYC': '#B85A00', 'BLINK_Kansas': '#1565C0',
    'SLEMM': '#222222', 'GCTA_MLMA': '#AA4499',
}
METHOD_MARKERS = dict(zip(METHOD_ORDER, ['o', 's', 'D', '^', 'v', 'P', '<', '>', 'X', 'p']))
METHOD_LINESTYLES = {key: ('--' if key in ['REGENIE_WO_LOCO', 'FarmCPU_Kansas',
                                          'BLINK_Kansas', 'GCTA_MLMA'] else '-')
                     for key in METHOD_ORDER}
METHOD_SHORT = dict(zip(METHOD_ORDER, ['BOLT', 'REG-L', 'REG-W', 'fastGWA',
                                      'FC-N', 'FC-K', 'BL-N', 'BL-K', 'SLEMM', 'GCTA']))
COLORS, LABELS, MARKERS, LINESTYLES = METHOD_COLORS, METHOD_LABELS, METHOD_MARKERS, METHOD_LINESTYLES


def white_contrast(color):
    rgb = [int(color[index:index+2], 16) / 255.0 for index in (1, 3, 5)]
    linear = [value / 12.92 if value <= .04045 else ((value + .055) / 1.055) ** 2.4
              for value in rgb]
    luminance = sum(weight * value for weight, value in zip((.2126, .7152, .0722), linear))
    return 1.05 / (luminance + .05)


def contrast_safe_ink(color):
    """Retain the plotted hue while making colored legend text legible on white."""
    rgb = [int(color[index:index+2], 16) for index in (1, 3, 5)]
    current = color
    while white_contrast(current) < 4.5:
        rgb = [max(0, int(round(value * .95))) for value in rgb]
        current = '#' + ''.join('{:02X}'.format(value) for value in rgb)
    return current


METHOD_LABEL_COLORS = {key: contrast_safe_ink(color) for key, color in METHOD_COLORS.items()}
REGION_COLORS = {'effect': '#D55E00', 'null': '#D9D9D9', 'central_null': '#8F8F8F'}
PIPELINE_COLORS = {'SuSiE-adj': '#009E73', 'SuSiE-RSS': '#0072B2'}
PIPELINE_MARKERS = {'SuSiE-adj': 'o', 'SuSiE-RSS': 's'}
PIPELINE_LINESTYLES = {'SuSiE-adj': '-', 'SuSiE-RSS': '--'}


def configure(fontsize=8.5):
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': fontsize,
        'axes.labelsize': fontsize, 'axes.titlesize': fontsize,
        'xtick.labelsize': fontsize, 'ytick.labelsize': fontsize,
        'legend.fontsize': fontsize, 'legend.title_fontsize': fontsize,
        'axes.linewidth': 0.65, 'lines.linewidth': 1.0,
        'lines.markersize': 3.5, 'lines.markeredgewidth': 0.45,
        'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'grid.linewidth': 0.35, 'grid.alpha': 0.25,
        'text.color': '#222222', 'axes.labelcolor': '#222222',
        'xtick.color': '#222222', 'ytick.color': '#222222',
        'axes.edgecolor': '#444444', 'axes.spines.top': False,
        'axes.spines.right': False, 'legend.frameon': False,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'savefig.facecolor': 'white', 'savefig.dpi': 300,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'pdf.use14corefonts': False,
        'svg.fonttype': 'none', 'mathtext.fontset': 'dejavusans',
        'pdf.compression': 9,
    })


apply_style = configure


def panel_label(ax, label, x=-0.12, y=1.04, **kwargs):
    options = dict(transform=ax.transAxes, fontweight='bold', fontsize=10,
                   ha='left', va='bottom', clip_on=False)
    options.update(kwargs)
    return ax.text(x, y, label.lower(), **options)


def save_figure(fig, stem, width_mm=170, height_mm=None, out_root=None, **kwargs):
    """Export one composite page at exact physical width; caller sets layout."""
    destination = Path(out_root) if out_root else ROOT
    for directory in ('pdfs', 'png'):
        (destination / directory).mkdir(parents=True, exist_ok=True)
    if height_mm is None:
        height_mm = fig.get_figheight() * 25.4
    if height_mm > 190.001:
        raise ValueError('Graphic height exceeds 190 mm: {}'.format(height_mm))
    fig.set_size_inches(width_mm / 25.4, height_mm / 25.4, forward=True)
    pdf = destination / 'pdfs' / (str(stem) + '.pdf')
    png = destination / 'png' / (str(stem) + '.png')
    fig.savefig(str(pdf), format='pdf', bbox_inches=None,
                metadata={'Title': str(stem), 'Creator': 'Matplotlib; verified benchmark source data'}, **kwargs)
    fig.savefig(str(png), format='png', dpi=300, bbox_inches=None, **kwargs)
    return pdf, png


configure()
