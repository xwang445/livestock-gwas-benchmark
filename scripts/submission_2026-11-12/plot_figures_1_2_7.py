#!/usr/bin/env python3
"""Publication redraws of LD decay, genomic design, and retained theory bins.

Figure 7 uses the original retained aggregate predictions and observations.
The unavailable variant-level theory generator is not reconstructed here.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from figure_style import (ROOT, DATA_ROOT, COLORS, LABELS, MARKERS, REGION_COLORS,
                          configure, panel_label, save_figure)

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


SOURCES = DATA_ROOT
CAPTIONS = ROOT / 'captions'
AUDIT = ROOT / 'audit'

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def caption(number, title, text):
    assert len(title.split()) <= 15
    assert len(text.split()) <= 300
    (CAPTIONS / 'Figure_{:02d}_caption.txt'.format(number)).write_text(
        'Figure {}. {}\n\n{}\n'.format(number, title, text))


def plot_ld():
    binned_path = SOURCES / 'figure01_ld_10kb_bins.tsv'
    curve_path = SOURCES / 'figure01_shared_spline.tsv'
    bins = pd.read_csv(binned_path, sep="\t")
    curve = pd.read_csv(curve_path, sep="\t")
    assert len(bins) == 1000 and bins.bin_id.is_unique
    assert np.allclose(bins.bin_end - bins.bin_start, .01)
    assert np.allclose(bins.mid_mb, (bins.bin_start + bins.bin_end) / 2)
    assert abs(bins.iloc[0].mean_r2 - .592525238698875) < 1e-12
    assert abs(curve.iloc[0].mean_r2 - bins.iloc[0].mean_r2) < 1e-12
    assert np.isfinite(curve[['mid_mb', 'mean_r2']].to_numpy()).all()
    fig, axes = plt.subplots(1, 2, figsize=(170 / 25.4, 73 / 25.4), sharey=True)
    fig.subplots_adjust(left=.092, right=.970, bottom=.19, top=.90, wspace=.15)
    for ax, maximum, letter in zip(axes, [2, 10], ['a', 'b']):
        view = bins.loc[bins.mid_mb <= maximum]
        fitted = curve.loc[curve.mid_mb <= maximum]
        ax.scatter(view.mid_mb, view.mean_r2, color='#9E9E9E', s=5.5,
                   linewidths=0, alpha=.65, zorder=2)
        ax.plot(fitted.mid_mb, fitted.mean_r2, color='#222222', linewidth=1.25,
                zorder=3)
        ax.set_xlim(0, maximum)
        ax.set_ylim(0, .62)
        ax.set_xticks(np.linspace(0, maximum, 5))
        ax.set_yticks(np.arange(0, .61, .1))
        ax.set_xlabel('Distance (Mb)')
        ax.grid(axis='y')
        panel_label(ax, letter, x=-.16)
        handles = [Line2D([], [], color='#222222', lw=1.25, label='Shared fitted curve'),
                   Line2D([], [], ls='none', marker='o', color='#9E9E9E', ms=3,
                          label='10 kb bin means')]
        ax.legend(handles=handles, loc='upper right', handlelength=1.7,
                  borderaxespad=.25, labelspacing=.35)
    axes[0].set_ylabel(r'Mean pairwise $r^2$')
    save_figure(fig, 'Figure_01', height_mm=73, out_root=ROOT)
    plt.close(fig)
    caption(1, 'Linkage disequilibrium decay across the pig genome',
            'Mean pairwise r² was calculated in identical 10 kb distance bins from the '
            '135,358-SNP panel in all 31,301 animals. Each retained SNP pair receives '
            'equal weight within its bin. Gray points show bin means; the black line '
            'is the shared fitted curve. The spline was fitted once to the complete '
            '0–10 Mb set of bin means, weighting bins equally, and reused for both '
            'views: (a) 0–2 Mb and (b) 0–10 Mb. Bins are left-closed and right-open, '
            'with 10 Mb retained in the last bin. The first bin mean is 0.5925 at '
            'the 5 kb midpoint, rather than an estimate at zero distance.')
    return dict(figure=1, sources=[dict(path=p.name, sha256=sha256(p))
                                  for p in [binned_path, curve_path]],
                n_bins=len(bins), retained_snp_pairs=int(bins.n_pairs.sum()),
                first_bin_mean=float(bins.iloc[0].mean_r2),
                same_spline_in_both_views=True)


def bracket(ax, x1, x2, y, text, fontsize=8.5):
    ax.plot([x1, x1, x2, x2], [y-.08, y, y, y-.08], color='#444444', lw=.8)
    ax.text((x1+x2)/2, y+.045, text, ha='center', va='bottom', fontsize=fontsize)


def plot_design():
    fig, ax = plt.subplots(figsize=(170 / 25.4, 105 / 25.4))
    fig.subplots_adjust(left=.04, right=.98, bottom=.045, top=.96)
    ax.set_xlim(-.22, 8.15)
    ax.set_ylim(-3.05, 2.20)
    ax.axis('off')
    ax.text(-.20, 1.95, 'a', fontsize=10, fontweight='bold', va='top')
    top, height = 1.02, .62
    pieces = [(0, 2, 'effect', 'Polygenic effect\nwindow (2 Mb)'),
              (2, 1, 'null', 'Null flank\n(1 Mb)'),
              (3, 2, 'central_null', 'Central-null\nregion (2 Mb)'),
              (5, 1, 'null', 'Null flank\n(1 Mb)'),
              (6, 2, 'effect', 'Polygenic effect\nwindow (2 Mb)')]
    for x, width, region, label in pieces:
        ax.add_patch(Rectangle((x, top), width, height,
                               facecolor=REGION_COLORS[region], edgecolor='#444444', lw=.7,
                               hatch='///' if region == 'central_null' else None))
        ax.text(x+width/2, top+height/2, label, ha='center', va='center',
                fontsize=8.5, color='#111111',
                bbox=dict(facecolor=REGION_COLORS[region], edgecolor='none', pad=.6))
    bracket(ax, 2, 6, 1.92, 'Null window (4 Mb)')
    for x in [0, 2, 3, 5, 6, 8]:
        ax.text(x, .91, '8 Mb' if x == 8 else '{}'.format(x),
                ha='right' if x == 8 else 'center', va='top', fontsize=8.5)
    # Exact 500 kb windows: the central-null region comprises four such windows.
    for i in range(16):
        x = i * .5
        region = 'effect' if x < 2 or x >= 6 else ('central_null' if 3 <= x < 5 else 'null')
        ax.add_patch(Rectangle((x, .12), .5, .34, facecolor=REGION_COLORS[region],
                               edgecolor='#555555', lw=.65,
                               hatch='///' if region == 'central_null' else None))
    ax.text(4, .58, 'Non-overlapping 500 kb locus windows', ha='center', va='bottom', fontsize=8.5)
    bracket(ax, 3, 5, -.20, 'Four central-null windows of 500 kb')
    ax.text(-.20, -.55, 'b', fontsize=10, fontweight='bold', va='top')
    ax.text(4, -.52, 'Fine-scale analysis: 80 consecutive 50 kb bins per null window',
            ha='center', va='top', fontsize=8.5)
    start, end = .35, 7.65
    cell_width = (end-start)/12
    indices = ['1', '2', '3', '…', '39', '40', '40', '39', '…', '3', '2', '1']
    for i, label in enumerate(indices):
        x = start+i*cell_width
        ax.add_patch(Rectangle((x, -1.62), cell_width, .43, facecolor='#F4F4F4',
                               edgecolor='#555555', lw=.65))
        ax.text(x+cell_width/2, -1.405, '…' if label == '…' else '50 kb',
                ha='center', va='center', fontsize=8.5)
        ax.text(x+cell_width/2, -1.76, label, ha='center', va='top', fontsize=8.5)
    ax.plot([4, 4], [-1.05, -2.32], color='#555555', ls='--', lw=.7)
    ax.plot([2, start], [.05, -1.14], ls='--', lw=.65, color='#777777')
    ax.plot([6, end], [.05, -1.14], ls='--', lw=.65, color='#777777')
    for x1, x2 in [(start, 3.96), (4.04, end)]:
        ax.annotate('', xy=(x1, -2.20), xytext=(x2, -2.20),
                    arrowprops=dict(arrowstyle='<->', color='#444444', lw=.8))
    ax.text((start+4)/2, -2.36, 'Index 1–40 from the left\neffect-window boundary',
            ha='center', va='top', fontsize=8.5)
    ax.text((4+end)/2, -2.36, 'Index 1–40 from the right\neffect-window boundary',
            ha='center', va='top', fontsize=8.5)
    save_figure(fig, 'Figure_02', height_mm=105, out_root=ROOT)
    plt.close(fig)
    caption(2, 'Genomic windows for effect and null-region analyses',
            '(a) The simulation alternates 2 Mb polygenic effect windows and 4 Mb '
            'null windows. Causal positions occur only in effect windows. Each null '
            'window contains a central 2 Mb region separated from either adjacent '
            'effect window by 1 Mb. Genome-wide detected loci are non-overlapping '
            '500 kb windows; each complete central-null region contains four such '
            'windows. The illustrated 8 Mb segment includes two adjacent effect '
            'windows. (b) For fine-scale spatial analyses, each 4 Mb null window '
            'is partitioned into 80 consecutive 50 kb bins. Bins are indexed '
            'symmetrically from 1 to 40 by their distance from the nearest '
            'effect-window boundary. Ellipses omit intervening bins; the expanded '
            'bin diagram is schematic. Central-null regions are nested within null '
            'regions, rather than forming an additional disjoint region class.')
    return dict(figure=2, source='Code-native diagram of verified analysis design',
                effect_window_mb=2, null_window_mb=4, central_null_mb=2,
                flank_mb=1, locus_window_kb=500, spatial_bin_kb=50,
                spatial_bins_per_null_window=80)


def plot_equation13():
    source = SOURCES / 'figure07_cached_theory_bin_summary.tsv'
    per_source = SOURCES / 'figure07_cached_observed_bins_by_phenotype.tsv'
    summary = pd.read_csv(source, sep='\t')
    per = pd.read_csv(per_source, sep='\t')
    methods = ['BOLT-LMM', 'REGENIE_LOCO']
    summary = summary.loc[(summary.region == 'Genome-wide') & summary.method.isin(methods)].copy()
    per = per.loc[(per.region == 'Genome-wide') & per.method.isin(methods)].copy()
    assert len(summary) == 40 and len(per) == 400
    check = per.groupby(['method', 'bin']).mean_observed_chisq.agg(['mean', 'std', 'size'])
    aligned = summary.set_index(['method', 'bin']).sort_index()
    assert np.allclose(aligned.mean_observed_chisq, check['mean'])
    assert np.allclose(aligned.sd_observed_chisq, check['std'])
    assert (check['size'] == 10).all()
    assert np.allclose(aligned.ci_upper - aligned.mean_observed_chisq,
                       2.2621571627409915*check['std']/np.sqrt(10), atol=1e-8)
    fig, axes = plt.subplots(1, 2, figsize=(170/25.4, 100/25.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.095, right=.985, top=.94, bottom=.19, wspace=.14)
    for ax, method, letter in zip(axes, methods, ['a', 'b']):
        q = summary.loc[summary.method == method].sort_values('bin')
        ax.plot([0, 70], [0, 70], color='#777777', ls='--', lw=.8, zorder=1)
        ax.errorbar(q.predicted_chisq, q.mean_observed_chisq,
                    yerr=[q.mean_observed_chisq-q.ci_lower, q.ci_upper-q.mean_observed_chisq],
                    fmt=MARKERS[method], ms=3.7, color=COLORS[method],
                    ecolor=COLORS[method], elinewidth=.8, capsize=1.7,
                    markeredgewidth=.5, zorder=3)
        ax.set_xlim(0, 70)
        ax.set_ylim(0, 70)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xticks([0, 20, 40, 60])
        ax.set_yticks([0, 20, 40, 60])
        ax.grid(axis='both')
        panel_label(ax, letter, x=-.16)
        handles = [Line2D([], [], color=COLORS[method], ls='none', marker=MARKERS[method],
                          label=LABELS[method], ms=3.7),
                   Line2D([], [], color='#777777', ls='--', lw=.8, label='Identity line')]
        ax.legend(handles=handles, loc='upper left', labelspacing=.4, handlelength=1.5)
    axes[0].set_ylabel(r'Observed mean $\chi^2$')
    fig.text(.54, .055, r'Predicted $\chi^2 = 1 + \mathrm{NCP}$', ha='center', va='bottom', fontsize=8.5)
    fig.canvas.draw()
    dimensions = []
    for ax in axes:
        box = ax.get_window_extent()
        width_mm, height_mm = box.width / fig.dpi * 25.4, box.height / fig.dpi * 25.4
        assert abs(width_mm - height_mm) < 1e-6
        dimensions.append(dict(x_axis_mm=width_mm, y_axis_mm=height_mm))
    save_figure(fig, 'Figure_07', height_mm=100, out_root=ROOT)
    plt.close(fig)
    caption(7, 'Predicted and observed association statistics for LOCO methods',
            'The 135,358 SNPs were grouped into 20 equal-count bins using the '
            'cached theoretical non-centrality parameter (NCP) from Equation 13. '
            'Each point compares the mean predicted statistic, 1 + NCP, with the '
            'observed mean association χ² for (a) BOLT-LMM and (b) REGENIE_LOCO. '
            'Observed values average the ten phenotype-specific bin means. '
            'Vertical error bars are 95% confidence intervals calculated as the '
            'mean ± t₉,₀.₉₇₅ × sample SD/√10. The dashed line is the identity line. '
            'Each panel has equal physical x- and y-axis lengths and common 0–70 limits. '
            'Predictions and bin assignments are retained from the original cached '
            'analysis; the figure is redrawn from the saved aggregate values.')
    return dict(figure=7, sources=[dict(path=p.name, sha256=sha256(p))
                                  for p in [source, per_source]],
                n_methods=2, n_bins_per_method=20, n_phenotypes=10,
                equal_axis_lengths_verified=True, axis_dimensions_mm=dimensions,
                aggregate_means_sd_and_ci_reproduced=True,
                ncp_predictions='Retained original cached values; upstream variant-level generator not rerun')


def main():
    global ROOT, SOURCES, CAPTIONS, AUDIT
    parser = argparse.ArgumentParser()
    parser.add_argument('--figure', type=int, choices=[1, 2, 7], help='Redraw one figure only')
    parser.add_argument('--source', type=Path, default=DATA_ROOT)
    parser.add_argument('--output', type=Path, default=ROOT)
    args = parser.parse_args()
    ROOT, SOURCES = args.output.resolve(), args.source.resolve()
    CAPTIONS, AUDIT = ROOT / 'captions', ROOT / 'audit'
    configure(8.5)
    for path in [SOURCES, CAPTIONS, AUDIT]:
        path.mkdir(parents=True, exist_ok=True)
    functions = {1: plot_ld, 2: plot_design, 7: plot_equation13}
    if args.figure is None:
        records = [functions[number]() for number in [1, 2, 7]]
    else:
        existing = AUDIT / 'figures_01_02_07_sources.json'
        previous = json.loads(existing.read_text())['figures'] if existing.exists() else []
        records = [record for record in previous if record['figure'] != args.figure]
        records.append(functions[args.figure]())
        records.sort(key=lambda record: record['figure'])
    payload = dict(figures=records, source_sha256=sha256(__file__),
                   output_width_mm=170, font_min_pt=8.5, min_linewidth_pt=.35,
                   titles_and_captions_separate=True)
    (AUDIT/'figures_01_02_07_sources.json').write_text(json.dumps(payload, indent=2)+'\n')
    print('Exported Figures 1, 2 and 7 as 170 mm embedded-font PDFs and 300 dpi PNGs.')


if __name__ == '__main__':
    main()
