#!/usr/bin/env python3
"""Export Figures 3, 5, 6 and ten full-size overlap supplements.

Sources are the current, inclusive-threshold canonical association tables.
The historical 26-locus SLEMM configuration is deliberately excluded.
Run: python3 scripts/plot_figures_3_5_6.py
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gwas_submission_matplotlib_20261112')
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter, ScalarFormatter

from figure_style import (
    ROOT, DATA_ROOT, METHOD_ORDER, LABELS, COLORS, MARKERS, LINESTYLES,
    configure, panel_label, save_figure, plt,
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def axes_mm(fig, left, bottom, width, height):
    w, h = fig.get_size_inches() * 25.4
    return fig.add_axes([left / w, bottom / h, width / w, height / h])


def text_color(rgba):
    rgb = np.asarray(rgba[:3])
    linear = np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
    luminance = linear @ np.array([.2126, .7152, .0722])
    return '#222222' if luminance > .179 else 'white'


def check_and_save(fig, stem, out, height, audit):
    # Keep the physical page fixed; tight cropping would change type sizes.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside, small = [], []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        if text.get_fontsize() < 8 - 1e-6:
            small.append([text.get_text(), text.get_fontsize()])
        extent = text.get_window_extent(renderer)
        if extent.width and extent.height and (
            extent.x0 < -.5 or extent.y0 < -.5 or
            extent.x1 > fig.bbox.width + .5 or extent.y1 > fig.bbox.height + .5
        ):
            outside.append(text.get_text())
    require(not outside, '{}: text outside page: {}'.format(stem, outside))
    require(not small, '{}: font below 8 pt: {}'.format(stem, small))
    save_figure(fig, stem, width_mm=170, height_mm=height, out_root=out)
    audit.append(dict(figure=stem, width_mm=170, height_mm=height,
                      minimum_font_pt=min(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text)
                                          if t.get_visible() and t.get_text().strip()),
                      text_outside_page=outside, vector_pdf_fonttype=42))
    plt.close(fig)


def caption(out, stem, text):
    target = out / 'captions' / (stem + '_caption.txt')
    target.parent.mkdir(parents=True, exist_ok=True)
    require(len(text.split()) <= 300, 'Caption exceeds 300 words: ' + stem)
    target.write_text(text.rstrip() + '\n')


def colorbar_strokes(colorbar):
    # Matplotlib 3.3 creates a white background Polygon with a 0.01-pt
    # edge. It is occluded by the color mesh, but normalize it as well so
    # every emitted PDF stroke exceeds the journal's 0.25-pt minimum.
    for name in ('patch', '_patch'):
        patch = getattr(colorbar, name, None)
        if patch is not None:
            patch.set_linewidth(.35)
    colorbar.ax.patch.set_linewidth(.35)


def source_tables(source, out):
    counts = pd.read_csv(source / 'Figure_03_counts.tsv', sep='\t')
    overlap = pd.read_csv(source / 'Figure_05_jaccard_by_phenotype.tsv', sep='\t')
    summary = pd.read_csv(source / 'Figure_06_ab_summary.tsv', sep='\t').set_index('configuration').reindex(METHOD_ORDER)
    means = pd.read_csv(source / 'Figure_06_cd_summary.tsv', sep='\t')
    require(len(counts) == 100 and not counts.duplicated(['configuration', 'phenotype']).any(),
            'Expected ten configurations by ten phenotypes.')
    totals = counts.groupby('configuration').loci.sum()
    require((totals['SLEMM'], totals['GCTA_MLMA']) == (16, 15), 'Sparse-method locus counts changed.')
    require(len(overlap) == 1000 and not overlap.duplicated(['phenotype', 'configuration_a', 'configuration_b']).any(),
            'Expected 1,000 phenotype-specific overlap values.')
    require((overlap.intersection <= overlap.union).all(), 'Invalid overlap counts.')
    defined = overlap.n_a.gt(0) & overlap.n_b.gt(0)
    require(overlap.loc[~defined, 'jaccard'].isna().all(), 'Empty-set Jaccard must be undefined.')
    require(np.allclose(overlap.loc[defined, 'jaccard'],
                        overlap.loc[defined, 'intersection'] / overlap.loc[defined, 'union']), 'Jaccard mismatch.')
    require(len(means) == 400 and means.n_phenotypes.eq(10).all(), 'Incomplete spatial summaries.')
    return counts, overlap, summary, means


def draw_counts(counts, out, audit):
    n = counts.pivot(index='configuration', columns='phenotype', values='loci').reindex(METHOD_ORDER)
    s = counts.pivot(index='configuration', columns='phenotype', values='significant_snps').reindex(METHOD_ORDER)
    fig = plt.figure(figsize=(170 / 25.4, 140 / 25.4))
    ax = axes_mm(fig, 38, 30, 126, 105)
    cmap, norm = plt.get_cmap('plasma'), LogNorm(vmin=1, vmax=float(n.to_numpy().max() + 1))
    mesh = ax.pcolormesh(np.arange(11), np.arange(11), n.to_numpy() + 1,
                         cmap=cmap, norm=norm, edgecolors='white', linewidth=.4, rasterized=False)
    for i in range(10):
        for j in range(10):
            foreground = text_color(cmap(norm(n.iloc[i, j] + 1)))
            ax.text(j + .5, i + .34, str(int(n.iloc[i, j])), ha='center', va='center',
                    fontsize=8.5, fontweight='bold', color=foreground)
            ax.text(j + .5, i + .69, '({:d})'.format(int(s.iloc[i, j])), ha='center', va='center',
                    fontsize=8.5, color=foreground)
    ax.set(xticks=np.arange(10) + .5, yticks=np.arange(10) + .5,
           xlim=(0, 10), ylim=(10, 0), xlabel='Phenotype')
    ax.set_xticklabels(range(1, 11))
    ax.set_yticklabels([LABELS[m] for m in METHOD_ORDER])
    ax.tick_params(axis='y', length=0, pad=5)
    ax.tick_params(axis='x', length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Center the colorbar under the main heatmap (38 + 126/2 =101mm),
    # with its descriptive key immediately to the left.
    cax = axes_mm(fig, 59, 9, 84, 3)
    cb = fig.colorbar(mesh, cax=cax, orientation='horizontal', ticks=[1, 11, 101, 1001])
    colorbar_strokes(cb)
    cb.ax.set_xticklabels(['0', '10', '100', '1,000'])
    cb.ax.minorticks_off()
    cb.solids.set_rasterized(False)
    fig.text(55 / 170, 10.5 / 140, 'Number of Detected loci', ha='right', va='center')
    check_and_save(fig, 'Figure_03', out, 140, audit)
    caption(out, 'Figure_03',
            'Figure 3. Significant loci and SNPs across simulated phenotypes. '
            'Each cell shows the number of detected loci above the number of significant SNPs in parentheses. '
            'Loci are chromosome-specific, zero-anchored 500-kb windows containing at least one SNP with '
            'P ≤ 5 × 10⁻⁸; terminal fragments are retained. Color represents detected-locus count, using '
            'a logarithmic scale for count plus one. All ten configurations use the same 20,000 animals, '
            '135,358 markers and ten phenotype replicates. SLEMM uses the updated benchmark configuration. '
            'REGENIE_LOCO and REGENIE_woLOCO denote analyses with and without '
            'leave-one-chromosome-out prediction, respectively. NYC and Kansas identify the corresponding '
            'FarmCPU and BLINK implementations. The source table contains all 100 count pairs.')


def overlap_matrix(overlap, phenotype):
    return overlap[overlap.phenotype.eq(phenotype)].pivot(
        index='configuration_a', columns='configuration_b', values='jaccard'
    ).reindex(index=METHOD_ORDER, columns=METHOD_ORDER).to_numpy()


def heatmap(ax, values, cmap, annotate=False):
    mesh = ax.pcolormesh(np.arange(11), np.arange(11), np.ma.masked_invalid(values),
                        cmap=cmap, vmin=0, vmax=1, edgecolors='white', linewidth=.35,
                        shading='flat', rasterized=False)
    ax.set(xticks=np.arange(10) + .5, yticks=np.arange(10) + .5, xlim=(0, 10), ylim=(10, 0))
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if annotate:
        for i in range(10):
            for j in range(10):
                value = values[i, j]
                foreground = '#222222' if pd.isna(value) else text_color(cmap(value))
                ax.text(j + .5, i + .5, 'NA' if pd.isna(value) else '{:.2f}'.format(value),
                        ha='center', va='center', fontsize=8.5, color=foreground)
    return mesh


def draw_overlap(overlap, out, audit):
    cmap = copy.copy(plt.get_cmap('Blues'))
    cmap.set_bad('#D9D9D9')
    fig = plt.figure(figsize=(170 / 25.4, 190 / 25.4))
    for phenotype in range(1, 11):
        row, col = divmod(phenotype - 1, 2)
        ax = axes_mm(fig, 35 + col * 67, 157.3 - row * 32.5, 61, 29)
        mesh = heatmap(ax, overlap_matrix(overlap, phenotype), cmap)
        ax.set_xticklabels([LABELS[m] for m in METHOD_ORDER] if row == 4 else [''] * 10,
                           rotation=90, ha='center', va='top', fontsize=8)
        ax.set_yticklabels([LABELS[m] for m in METHOD_ORDER] if col == 0 else [''] * 10,
                           fontsize=8)
        ax.tick_params(axis='x', pad=2)
        ax.tick_params(axis='y', pad=4)
        ax.text(0, 1.005, 'Phenotype {}'.format(phenotype),
                transform=ax.transAxes, fontsize=8.5, ha='left', va='bottom')
    fig.text(2 / 170, 108 / 190, 'Method', rotation=90, ha='center', va='center', fontsize=8.5)
    fig.text(.58, 5 / 190, 'Method', ha='center', va='center', fontsize=8.5)
    cax = axes_mm(fig, 5, 4, 25, 2.3)
    cb = fig.colorbar(mesh, cax=cax, orientation='horizontal', ticks=[0, .5, 1])
    colorbar_strokes(cb)
    cb.solids.set_rasterized(False)
    cb.ax.tick_params(labelsize=8, pad=1, length=1.5)
    fig.text(17.5 / 170, 10.5 / 190, 'Jaccard', ha='center', va='center', fontsize=8)
    fig.legend(handles=[Patch(facecolor='#D9D9D9', edgecolor='#777777', linewidth=.4,
                              label='NA')],
               loc='lower left', bbox_to_anchor=(5 / 170, 16 / 190), fontsize=8,
               borderaxespad=0, handlelength=1.0)
    check_and_save(fig, 'Figure_05', out, 190, audit)
    caption(out, 'Figure_05',
            'Figure 5. Phenotype-specific overlap of detected loci. '
            'Each panel corresponds to the indicated phenotype. Method names on the left apply to '
            'all panels, and the shared vertical labels beneath the bottom panels give the column '
            'methods. Cell color represents the Jaccard '
            'index, the number of shared detected 500-kb windows divided by their union, using '
            'P ≤ 5 × 10⁻⁸. Grey cells are undefined because either configuration has no detected loci; '
            'this convention also applies to empty-set diagonal cells. SLEMM and GCTA-MLMA have no '
            'detected loci in phenotypes 2, 6, 8, 9 and 10. Exact cell values and method names are '
            'provided in the source table and in full-size annotated matrices in Figures S11–S20.')
    supplement = out / 'supplement'
    for phenotype in range(1, 11):
        stem = 'Figure_S{:02d}'.format(10 + phenotype)
        fig = plt.figure(figsize=(170 / 25.4, 150 / 25.4))
        ax = axes_mm(fig, 38, 44, 126, 100)
        mesh = heatmap(ax, overlap_matrix(overlap, phenotype), cmap, annotate=True)
        ax.set_yticklabels([LABELS[m] for m in METHOD_ORDER])
        ax.set_xticklabels([LABELS[m] for m in METHOD_ORDER], rotation=50, ha='right')
        cax = axes_mm(fig, 81, 7, 76, 3)
        cb = fig.colorbar(mesh, cax=cax, orientation='horizontal', ticks=[0, .25, .5, .75, 1])
        colorbar_strokes(cb)
        cb.solids.set_rasterized(False)
        fig.text(75 / 170, 8.5 / 150, 'Jaccard', ha='right', va='center')
        check_and_save(fig, stem, supplement, 150, audit)
        caption(supplement, stem,
                'Figure S{}. Detected-locus overlap for phenotype {}. '
                'Cells show the Jaccard index between detected-locus sets for the indicated methods; '
                'loci are fixed, chromosome-specific 500-kb windows containing a SNP with P ≤ 5 × 10⁻⁸. '
                'Indices are rounded to two decimal places; 0.00 can therefore denote a nonzero value '
                'below 0.005. Grey cells marked NA are undefined because either locus set is empty. '
                'The calculation uses the current SLEMM configuration and all ten completed GCTA-MLMA '
                'phenotypes. Full-precision values are supplied in Figure_05_jaccard_by_phenotype.tsv.'
                .format(10 + phenotype, phenotype))


def draw_spillover(summary, means, out, audit):
    fig = plt.figure(figsize=(170 / 25.4, 160 / 25.4))
    a = axes_mm(fig, 39, 98, 43, 55)
    b = axes_mm(fig, 39, 29, 43, 50)
    c = axes_mm(fig, 108, 98, 58, 55)
    d = axes_mm(fig, 108, 29, 58, 50)
    y = np.arange(10)
    for ax in (a, b):
        ax.set_yticks(y)
        ax.set_yticklabels([LABELS[m] for m in METHOD_ORDER])
        ax.set_ylim(9.7, -.7)
        ax.grid(axis='x', color='#D9D9D9', linewidth=.35)
        ax.set_axisbelow(True)
        ax.tick_params(axis='y', length=0, pad=4)
    a.barh(y, summary.spillover_rate_mean * 100, xerr=summary.spillover_rate_sd * 100,
           color=[COLORS[m] for m in METHOD_ORDER], height=.64, edgecolor='white', linewidth=.4,
           error_kw=dict(elinewidth=.65, capsize=1.5, capthick=.65, ecolor='#333333'))
    for i, row in enumerate(summary.itertuples()):
        x = 100 * (row.spillover_rate_mean + row.spillover_rate_sd) + 1
        a.text(x, i, '{:.1f}%'.format(100 * row.spillover_rate_mean), va='center', fontsize=8.5)
    a.set_xlim(0, 45)
    a.set_xticks([0, 15, 30, 45])
    a.set_xlabel('Central-null windows detected (%)', labelpad=4)
    b.axvline(1, color='#555555', linewidth=.7, linestyle='--', zorder=1)
    for i, method in enumerate(METHOD_ORDER):
        b.errorbar(summary.loc[method, 'historical_lambda_mean'], i,
                   xerr=summary.loc[method, 'historical_lambda_sd'], color=COLORS[method],
                   marker=MARKERS[method], markersize=3.8, linestyle='none',
                   elinewidth=.8, capsize=1.5, capthick=.65)
    b.set_xlim(0, 10.5)
    b.set_xticks([0, 1, 5, 10])
    b.set_xlabel(r'Central-null $\lambda_{GC}$', labelpad=4)
    b.annotate(r'$\lambda_{GC}=1$', xy=(1, .04), xycoords=('data', 'axes fraction'),
               xytext=(4, 0), textcoords='offset points', ha='left', va='bottom', fontsize=8)
    for ax, methods in [(c, METHOD_ORDER), (d, METHOD_ORDER[2:])]:
        for method in methods:
            z = means[means.configuration.eq(method)].sort_values('distance_index')
            marker_every = (2, 5) if method == 'GCTA_MLMA' else (0, 5)
            ax.plot(z.distance_midpoint_mb, z.observed_expected_mean,
                    color=COLORS[method], marker=MARKERS[method], linestyle=LINESTYLES[method],
                    linewidth=1.0, markersize=3.2, markevery=marker_every,
                    markerfacecolor='white' if method == 'GCTA_MLMA' else COLORS[method],
                    markeredgewidth=.5, label=LABELS[method])
        ax.axhline(1, color='#777777', linestyle=':', linewidth=.65, zorder=0)
        ax.set_xlim(0, 2)
        ax.set_xticks([0, 1, 2])
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Distance from effect boundary (Mb)', labelpad=4)
        ax.set_ylabel('Observed / expected SNPs', labelpad=3)
        ax.grid(color='#D9D9D9', linewidth=.35)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((0, 0))
        ax.yaxis.set_major_formatter(formatter)
    for label, ax in [('a', a), ('b', b), ('c', c), ('d', d)]:
        panel_label(ax, label, x=-.07, y=1.02)
    handles = [Line2D([], [], color=COLORS[m], marker=MARKERS[m], linestyle=LINESTYLES[m],
                      linewidth=1.0, markersize=3.5, label=LABELS[m],
                      markerfacecolor='white' if m == 'GCTA_MLMA' else COLORS[m])
               for m in METHOD_ORDER]
    # Two compact rows, arranged in the shared method order when read row-wise.
    indices = [0, 5, 1, 6, 2, 7, 3, 8, 4, 9]
    fig.legend(handles=[handles[i] for i in indices], loc='lower center',
               bbox_to_anchor=(.52, .003), ncol=5, columnspacing=.9,
               handlelength=1.3, handletextpad=.4, labelspacing=.5, borderaxespad=0,
               fontsize=8.5)
    check_and_save(fig, 'Figure_06', out, 160, audit)
    caption(out, 'Figure_06',
            'Figure 6. Central-null spillover and spatial signal decay. '
            '(a) Fraction of 1,480 central-null 500-kb windows containing at least one SNP with '
            'P ≤ 5 × 10⁻⁸, shown as mean ± sample SD across ten phenotypes; labels give mean percentages. '
            '(b) Central-null genomic inflation, shown as mean ± sample SD across phenotypes. The historical '
            'Figure 6b estimator is retained: median finite one-degree-of-freedom association χ² divided '
            'by 0.4549364, using (β/SE)² when coefficients and standard errors are available. The dashed '
            'line marks λGC = 1. (c) Observed-to-expected significant-SNP counts in 50-kb null windows, '
            'with expected count equal to the tested-SNP count multiplied by 5 × 10⁻⁸. Ratios are averaged '
            'equally over nonempty windows at each mirrored distance index, then equally across ten '
            'phenotypes. Distances are bin midpoints measured from the nearer simulated effect boundary. '
            'The historical closed-endpoint convention is retained; source tables document boundary cases. '
            '(d) The same calculation excluding BOLT-LMM and REGENIE_LOCO, using its own vertical scale. '
            'Dotted reference lines in c and d indicate observed/expected = 1. Symbols and line styles '
            'supplement the shared method colors.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DATA_ROOT)
    parser.add_argument('--output', type=Path, default=ROOT)
    parser.add_argument('--figure', type=int, choices=[3, 5, 6], help='Redraw one main figure only')
    args = parser.parse_args()
    configure(8.5)
    out = args.output.resolve()
    for folder in ['source_data', 'captions', 'pdfs', 'png']:
        (out / folder).mkdir(parents=True, exist_ok=True)
    counts, overlap, summary, means = source_tables(args.source, out)
    audit = []
    if args.figure is not None:
        previous = out / 'source_data/Figures_03_05_06_provenance.json'
        previous_records = json.loads(previous.read_text())['figures'] if previous.exists() else []
        excluded = {'Figure_{:02d}'.format(args.figure)}
        if args.figure == 5:
            excluded.update('Figure_S{:02d}'.format(i) for i in range(11, 21))
        audit = [row for row in previous_records if row['figure'] not in excluded]
    if args.figure in (None, 3):
        draw_counts(counts, out, audit)
    if args.figure in (None, 5):
        draw_overlap(overlap, out, audit)
    if args.figure in (None, 6):
        draw_spillover(summary, means, out, audit)
    audit.sort(key=lambda row: row['figure'])
    sources = ['Figure_03_counts.tsv', 'Figure_05_jaccard_by_phenotype.tsv',
               'Figure_06_ab_summary.tsv', 'Figure_06_cd_summary.tsv']
    manifest = dict(status='complete', figures=audit, significance='P <= 5e-8',
                    canonical_slemm_loci=16, canonical_gcta_loci=15,
                    historical_SLEMM_previous_excluded=True,
                    fig6b_column='lambda_gc_central_null_historical',
                    fig6cd_column='observed_expected_equal_nonempty_bins',
                    source_sha256={name: digest(args.source / name) for name in sources},
                    script_sha256=digest(Path(__file__)))
    (out / 'source_data/Figures_03_05_06_provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(status='complete', n_exported=len(audit),
                          figures=[row['figure'] for row in audit]), indent=2))


if __name__ == '__main__':
    main()
