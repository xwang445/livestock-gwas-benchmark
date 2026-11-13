#!/usr/bin/env python3
"""Figure 4 and S1--S10 from all 100 current association outputs.

Preserves all ten phenotypes, with methods separated into 2/4/4 columns.
Sources are never modified. Source point tables contain variant coordinates
and association statistics, with no individual records.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/manhattan_submission_20261112')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
OUT = REPOSITORY / 'figures/submission_2026-11-12'
DATA_ROOT = REPOSITORY / 'derived_data/submission_2026-11-12'
BENCH = None
GCTA = None
ALPHA = 5e-8


def significant(values, is_log10=False):
    values = np.asarray(values)
    return values >= -np.log10(ALPHA) if is_log10 else values <= ALPHA


def source_config(bench, method, ph, gcta):
    """Explicit current sources; independent of any code outside this package."""
    if method == 'BOLT-LMM':
        return bench/('BOLTLMM/bolt_pheno%d.stats'%ph), r'\s+', ['CHR','BP','SNP','P_BOLT_LMM_INF'], False
    if method.startswith('REGENIE'):
        directory = 'REGENIE' if method == 'REGENIE_LOCO' else 'REGENIE_WO_LOCO'
        stem = 'regenie_pheno' if method == 'REGENIE_LOCO' else 'regenie_pheno_wo_loco'
        return bench/('%s/%s%d_pheno%d.regenie'%(directory,stem,ph,ph)), r'\s+', ['CHROM','GENPOS','ID','LOG10P'], True
    if method == 'fastGWA':
        return bench/('fastGWA/spg.even.100k.pheno%d.fastgwa.fastGWA'%ph), r'\s+', ['CHR','POS','SNP','P'], False
    if method == 'SLEMM':
        return bench/('SLEMM_qf90_window3000/spg.even.100k.pheno%d.allchr.gwas.txt'%ph), r'\s+', ['#CHROM','POS','ID','P'], False
    if method == 'GCTA_MLMA':
        return gcta/('spg.even.100k.pheno%d.mlma'%ph), r'\s+', ['Chr','bp','SNP','p'], False
    model, city = method.split('_')
    p = bench/('GAPIT/GAPIT.Association.GWAS_Results.%s.pheno%d(%s).csv'%(model,ph,city))
    return p, ',', ['Chr','Pos','SNP','P.value'], False

METHODS = ['BOLT-LMM', 'REGENIE_LOCO', 'REGENIE_WO_LOCO', 'fastGWA',
           'FarmCPU_NYC', 'FarmCPU_Kansas', 'BLINK_NYC', 'BLINK_Kansas',
           'SLEMM', 'GCTA_MLMA']
LABELS = dict(zip(METHODS, ['BOLT-LMM', 'REGENIE (LOCO)', 'REGENIE (no LOCO)',
    'fastGWA', 'FarmCPU (NYC)', 'FarmCPU (Kansas)', 'BLINK (NYC)',
    'BLINK (Kansas)', 'SLEMM', 'GCTA MLMA']))
COLORS = dict(zip(METHODS, ['#0072B2', '#D55E00', '#E69F00', '#009E73',
    '#D81B60', '#004D40', '#B85A00', '#1565C0', '#222222', '#AA4499']))
LABEL_COLORS = dict(COLORS)
LABEL_COLORS.update({'REGENIE_LOCO':'#9C4500','REGENIE_WO_LOCO':'#8C5E00',
    'fastGWA':'#006B4E'})
MARKERS = dict(zip(METHODS, ['o', '^', 's', 'D', '+', 'x', 'v', '*', 'o', 'p']))
THRESHOLD = -np.log10(ALPHA)


def configure_style():
    """Use the package defaults when available, then enforce export limits."""
    p = HERE / 'figure_style.py'
    if p.exists():
        spec = importlib.util.spec_from_file_location('submission_style', str(p))
        style = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(style)
        for name in ['apply_style', 'set_style', 'configure']:
            fn = getattr(style, name, None)
            if callable(fn):
                try:
                    fn()
                except TypeError:
                    pass
                break
        for name, target in [('COLORS', COLORS), ('METHOD_COLORS', COLORS),
                             ('LABELS', LABELS), ('METHOD_LABELS', LABELS),
                             ('MARKERS', MARKERS), ('METHOD_MARKERS', MARKERS),
                             ('METHOD_LABEL_COLORS', LABEL_COLORS)]:
            val = getattr(style, name, None)
            if isinstance(val, dict):
                for key in METHODS:
                    if key in val:
                        target[key] = val[key]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5,
        'axes.labelsize': 8.5, 'axes.titlesize': 8.5, 'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5, 'legend.fontsize': 8.5, 'pdf.fonttype': 42,
        'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': .5,
        'xtick.major.width': .5, 'ytick.major.width': .5,
        'xtick.major.size': 2, 'ytick.major.size': 2,
        'savefig.dpi': 300, 'figure.dpi': 150})


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def white_contrast(color):
    rgb=[int(color[i:i+2],16)/255. for i in (1,3,5)]
    linear=[v/12.92 if v<=.04045 else ((v+.055)/1.055)**2.4 for v in rgb]
    luminance=sum(w*v for w,v in zip((.2126,.7152,.0722),linear))
    return 1.05/(luminance+.05)


def load(method, phenotype):
    p, sep, cols, islog = source_config(BENCH, method, phenotype,
        GCTA)
    d = pd.read_csv(p, sep=sep, usecols=cols[:4])
    d = d[cols[:4]].copy()
    d.columns = ['chr', 'position', 'SNP', 'native_value']
    d['chr'] = pd.to_numeric(d.chr.astype(str).str.replace(r'^chr', '', regex=True))
    v = d.native_value.to_numpy()
    assert len(d) == 135358 and np.isfinite(v).all()
    assert d.chr.between(1, 18).all() and not d[['chr', 'position']].duplicated().any()
    assert (v >= 0).all() and (islog or (v <= 1).all())
    d['significant'] = significant(v, is_log10=islog)
    d['minus_log10p'] = v if islog else -np.log10(np.maximum(v, 1e-300))
    d = d.sort_values(['chr', 'position']).reset_index(drop=True)
    return d, p, islog


def clean_axis(ax):
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    for side in ['bottom', 'left']:
        ax.spines[side].set_color('#555555')
    ax.tick_params(pad=1.5)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_axisbelow(True)
    ax.grid(axis='y', color='#E8E8E8', linewidth=.35)


def check_and_export(fig, pdf, png, audit, high_resolution_png=None):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    outside = []
    sizes = []
    for obj in fig.findobj(match=matplotlib.text.Text):
        if not obj.get_visible() or not obj.get_text():
            continue
        box = obj.get_window_extent(renderer)
        sizes.append(float(obj.get_fontsize()))
        if box.x0 < -1 or box.y0 < -1 or box.x1 > width + 1 or box.y1 > height + 1:
            outside.append(obj.get_text())
    assert not outside, 'Text outside page: ' + repr(outside)
    assert min(sizes) >= 8
    pdf.parent.mkdir(parents=True, exist_ok=True)
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, dpi=300, facecolor='white', metadata={'Creator': 'Matplotlib; current benchmark outputs'})
    fig.savefig(png, dpi=300, facecolor='white')
    if high_resolution_png is not None:
        high_resolution_png.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(high_resolution_png,dpi=600,facecolor='white')
        assert high_resolution_png.stat().st_size < 10000000
    assert pdf.stat().st_size < 10000000 and png.stat().st_size < 10000000
    audit.append(dict(pdf=str(pdf.relative_to(OUT)), png=str(png.relative_to(OUT)),
        width_mm=round(fig.get_figwidth()*25.4, 6), height_mm=round(fig.get_figheight()*25.4, 6),
        min_font_pt=min(sizes), text_within_page=True, pdf_bytes=pdf.stat().st_size,
        png_bytes=png.stat().st_size, points_rasterized_at_dpi=300,
        axes_and_text_vector=True, minimum_nonzero_linewidth_pt=.35))
    plt.close(fig)


def main_figure(data, offsets, centers, xmax, audit):
    """One phenotype per row; separate 2/4/4-method columns reduce overlap."""
    fig = plt.figure(figsize=(170/25.4, 190/25.4))
    grid = fig.add_gridspec(10, 3, left=.20, right=.995, bottom=.18,
        top=.935, wspace=.17, hspace=.32)
    groups = [METHODS[:2], ['FarmCPU_NYC','FarmCPU_Kansas','BLINK_NYC','BLINK_Kansas'],
        ['fastGWA','REGENIE_WO_LOCO','SLEMM','GCTA_MLMA']]
    global_upper = float(np.ceil(max(d.minus_log10p.max() for d in data.values())/10)*10)
    upper = [global_upper]*3
    headings = ['LOCO methods', 'FarmCPU and BLINK', 'Other non-LOCO\nmethods']
    ticks = [0,5,11,17]
    axes = []
    for ph in range(1,11):
        row = []
        for col,group in enumerate(groups):
            ax = fig.add_subplot(grid[ph-1,col])
            row.append(ax)
            for method in group:
                d = data[(method,ph)]
                d = d[d.significant]
                options = dict(s=2.1, marker=MARKERS[method], color=COLORS[method],
                    linewidths=0, alpha=.8, rasterized=True)
                if method in groups[1]:
                    options['alpha'] = 1.0
                if method == 'SLEMM':
                    options.update(s=10, facecolors='none', edgecolors=COLORS[method],
                        linewidths=.6,zorder=6,alpha=1)
                elif method == 'GCTA_MLMA':
                    options.update(s=11,facecolors='none',edgecolors=COLORS[method],
                        linewidths=.6,zorder=7,alpha=1)
                ax.scatter((d.position+d.chr.map(offsets))/1e6,d.minus_log10p,**options)
            clean_axis(ax)
            ax.set(xlim=(0,xmax),ylim=(0,upper[col]),yticks=[0,upper[col]])
            ax.axhline(THRESHOLD,color='#777777',linestyle=(0,(3,2)),linewidth=.5,zorder=2)
            ax.set_xticks(centers.iloc[ticks]/1e6)
            ax.set_xticklabels([str(int(centers.index[i])) for i in ticks])
            ax.get_xticklabels()[-1].set_ha('right')
            if ph < 10:
                ax.tick_params(labelbottom=False,bottom=False)
        position = row[0].get_position()
        fig.text(.15,(position.y0+position.y1)/2,'Phenotype {}'.format(ph),
            ha='right',va='center',fontsize=8.5,fontweight='bold')
        axes.append(row)
    method_legends, heading_text = [], []
    for col,group in enumerate(groups):
        position = axes[0][col].get_position()
        center = (position.x0+position.x1)/2
        heading_text.append(fig.text(center,.953,headings[col],ha='center',va='bottom',
            fontsize=9.5,fontweight='medium'))
        handles = [Line2D([],[],marker=MARKERS[m],color=COLORS[m],linestyle='none',
            markerfacecolor='none' if m in ['SLEMM','GCTA_MLMA'] else COLORS[m],
            markeredgewidth=.6,markersize=4.5,label=LABELS[m]) for m in group]
        # The three key blocks align with the corresponding method columns.
        method_legends.append(fig.legend(handles,[LABELS[m] for m in group],ncol=1,loc='lower center',
            bbox_to_anchor=(center,.03),frameon=False,handlelength=1.,fontsize=10,
            prop={'size':10,'weight':'medium'},handletextpad=.35,labelspacing=.18,borderaxespad=0))
        for text, method in zip(method_legends[-1].get_texts(),group):
            text.set_color(LABEL_COLORS[method])
    threshold_key = fig.legend([Line2D([],[],color='#777777',linestyle=(0,(3,2)),linewidth=.5)],
        ['P = 5 × 10⁻⁸'],loc='lower center',bbox_to_anchor=(.56,.001),
        frameon=False,handlelength=1.5,borderaxespad=0)
    fig.text(.006,.55,'−log₁₀(P)',rotation=90,va='center',fontsize=8.5)
    chromosome_label = fig.text(.60,.14,'Chromosome',ha='center',fontsize=8.5)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    legend_boxes=[legend.get_window_extent(renderer) for legend in method_legends]
    protected=[chromosome_label.get_window_extent(renderer),threshold_key.get_window_extent(renderer)]
    protected.extend(text.get_window_extent(renderer) for ax in axes[-1]
        for text in ax.get_xticklabels() if text.get_visible())
    assert all(not a.overlaps(b) for a in legend_boxes for b in protected), 'Method legend intersects chromosome or significance labels'
    assert all(not legend_boxes[i].overlaps(legend_boxes[j]) for i in range(3) for j in range(i+1,3)), 'Method legend blocks overlap'
    assert all(tuple(ax.get_ylim())==(0.,global_upper) for row in axes for ax in row)
    assert all(white_contrast(LABEL_COLORS[m])>=4.5 for m in METHODS)
    assert all(white_contrast(COLORS[m])>=4.5 for m in groups[1])
    highres=OUT/'high_resolution/Figure_04_600dpi.png'
    label_audit=dict(status='PASS',layout='ten phenotype rows by three method-group columns',
        width_mm=170,height_mm=190,method_legend_font_sizes_pt=[text.get_fontsize()
            for legend in method_legends for text in legend.get_texts()],
        group_heading_font_sizes_pt=[text.get_fontsize() for text in heading_text],
        legend_marker_size_pt=4.5,primary_png_dpi=300,additional_png_dpi=600,
        high_resolution_png=str(highres.relative_to(OUT)),
        method_groups=groups,group_ordinate_ranges=[[0,v] for v in upper],
        method_point_colors={m:COLORS[m] for m in METHODS},
        method_label_colors={m:LABEL_COLORS[m] for m in METHODS},
        method_label_white_contrast={m:white_contrast(LABEL_COLORS[m]) for m in METHODS},
        middle_group_point_white_contrast={m:white_contrast(COLORS[m]) for m in groups[1]},
        all_30_axes_same_y_range=True,
        middle_group_point_opacity=1.0,
        method_legend_overlaps_chromosome_or_significance_labels=False,
        method_legend_blocks_overlap=False,all_30_axes_and_data_preserved=True,
        supplements_S01_S10_unchanged=True)
    check_and_export(fig,OUT/'pdfs/Figure_04.pdf',OUT/'png/Figure_04.png',audit,
        high_resolution_png=highres)
    label_audit.update(pdf_sha256=sha(OUT/'pdfs/Figure_04.pdf'),
        primary_png_sha256=sha(OUT/'png/Figure_04.png'),high_resolution_png_sha256=sha(highres),
        pdf_bytes=(OUT/'pdfs/Figure_04.pdf').stat().st_size,
        primary_png_bytes=(OUT/'png/Figure_04.png').stat().st_size,
        high_resolution_png_bytes=highres.stat().st_size)
    (OUT/'audit/figure4_label_resolution.json').write_text(json.dumps(label_audit,indent=2)+'\n')
    (OUT/'audit/figure4_contrast_grouping.json').write_text(json.dumps(label_audit,indent=2)+'\n')
    audit[-1].update(layout='ten phenotype rows by three method-group columns',
        group_ordinate_ranges=[[0,v] for v in upper],legend='bottom, aligned to method columns',
        phenotype_panel_letters=False)


def supplement(data, ph, offsets, centers, xmax, audit):
    fig, axes = plt.subplots(5,2,figsize=(170/25.4,175/25.4))
    fig.subplots_adjust(left=.095,right=.99,top=.925,bottom=.105,wspace=.19,hspace=.42)
    tick_ix = [0,2,4,6,8,10,12,14,17]
    ranges={};headers=[];source_maxima={};large_marker_counts={}
    for i,(method,ax) in enumerate(zip(METHODS,axes.flat)):
        d = data[(method,ph)]
        assert len(d)==135358
        xx = (d.position+d.chr.map(offsets))/1e6
        for parity,color in [(1,'#666666'),(0,'#A0A0A0')]:
            s = d.chr.mod(2).eq(parity)
            ax.scatter(xx[s],d.minus_log10p[s],s=.52,linewidths=0,alpha=1,
                color=color,rasterized=True)
        # Highlight the same significant SNPs already present in the complete
        # background. No tested variant is filtered out of these panels.
        selected=d.significant
        ax.scatter(xx[selected],d.minus_log10p[selected],s=.9,linewidths=0,
            alpha=1,color=COLORS[method],rasterized=True,zorder=3)
        if method in ['SLEMM','GCTA_MLMA']:
            large_marker_counts[method]=int(selected.sum())
            ax.scatter(xx[selected],d.minus_log10p[selected],s=16,marker='o',
                linewidths=.35,edgecolors='white',facecolors=COLORS[method],
                alpha=1,rasterized=True,zorder=4)
        ymax = 110.
        source_maxima[method]=float(d.minus_log10p.max())
        assert d.minus_log10p.ge(0).all() and source_maxima[method]<ymax
        ranges[method]=[0,ymax]
        ax.set(xlim=(0,xmax),ylim=(0,ymax))
        ax.set_yticks([0,55,110])
        ax.axhline(THRESHOLD,color='#777777',linestyle=(0,(3,2)),linewidth=.55,zorder=2)
        clean_axis(ax)
        headers.append(ax.text(.015,1.025,LABELS[method],transform=ax.transAxes,
            ha='left',va='bottom',fontsize=9.5,fontweight='medium',
            color=LABEL_COLORS[method],clip_on=False))
        ax.set_xticks(centers.iloc[tick_ix]/1e6)
        ax.set_xticklabels([str(int(centers.index[j])) for j in tick_ix])
        ax.get_xticklabels()[-1].set_ha('right')
        if i < 8:
            ax.tick_params(labelbottom=False,bottom=False)
    phenotype_label=fig.text(.54,.975,'Phenotype {}'.format(ph),ha='center',va='center',
        fontsize=10,fontweight='bold')
    fig.text(.006,.53,'−log₁₀(P)',rotation=90,va='center',fontsize=8.5)
    fig.text(.54,.018,'Chromosome',ha='center',fontsize=8.5)
    # Inset key lies in the reserved lower band; no manuscript caption is embedded.
    fig.legend([Line2D([],[],color='#777777',linestyle=(0,(3,2)),linewidth=.55)],
        ['P = 5 × 10⁻⁸'],loc='lower center',bbox_to_anchor=(.54,.045),frameon=False)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    data_boxes=[ax.get_window_extent(renderer) for ax in axes.flat]
    header_boxes=[text.get_window_extent(renderer) for text in headers]
    assert all(not label.overlaps(area) for label in header_boxes for area in data_boxes), 'Method header overlaps a data axis'
    assert all(not label.overlaps(phenotype_label.get_window_extent(renderer)) for label in header_boxes), 'Method header overlaps phenotype label'
    assert all(tuple(ax.get_ylim())==(0.,110.) for ax in axes.flat)
    check_and_export(fig,OUT/('supplement/pdfs/Figure_S%02d.pdf'%ph),
        OUT/('supplement/png/Figure_S%02d.png'%ph),audit)
    audit[-1].update(phenotype=ph,condition_label='Phenotype {}'.format(ph),
        method_headers=[LABELS[m] for m in METHODS],panel_letters=False,
        unique_tested_SNPs_per_panel=135358,n_method_panels=10,
        background='alternating gray chromosomes; all tested SNPs',
        significant_point_colors={m:COLORS[m] for m in METHODS},
        significant_highlights_duplicate_existing_background_SNPs=True,
        method_header_colors={m:LABEL_COLORS[m] for m in METHODS},
        ordinate_ranges=ranges,ordinate_policy='shared 0-110 across all 100 S1-S10 panels',
        source_maxima=source_maxima,observed_values_above_display_range=0,
        significant_marker_areas_pt2={m:(16 if m in ['SLEMM','GCTA_MLMA'] else .9) for m in METHODS},
        original_significant_point_layer_preserved=True,
        enlarged_significant_marker_shape='circle',enlarged_marker_white_edge_pt=.35,
        enlarged_significant_marker_counts=large_marker_counts,
        enlarged_marker_fill_colors={m:COLORS[m] for m in large_marker_counts},
        method_headers_above_axes=True,method_headers_overlap_data=False,
        header_background_boxes=False,
        chromosome_coordinates='shared chromosome 1-18 offsets across all panels')


def write_captions(main_only=False,supplement_only=False):
    if not supplement_only:
        (OUT/'captions/Figure_04_caption.txt').write_text(
        'Figure 4. Significant association signals across ten simulated phenotypes. '
        'Rows show phenotypes 1–10. Columns separate BOLT-LMM and REGENIE_LOCO; '
        'FarmCPU (NYC), FarmCPU (Kansas), BLINK (NYC) and BLINK (Kansas); and '
        'fastGWA, REGENIE_woLOCO, current SLEMM and GCTA-MLMA, respectively. '
        'All significant SNPs (P ≤ 5 × 10⁻⁸) are shown without truncation. '
        'All 30 axes use the same ordinate range, 0–110. The colored bottom keys '
        'identify configurations. Dashed '
        'lines mark the significance threshold. Full-genome plots containing all '
        'tested SNPs are supplied in Figures S1–S10.\n')
    if main_only:
        return
    for ph in range(1,11):
        (OUT/('supplement/captions/Figure_S%02d_caption.txt'%ph)).write_text(
            'Figure S%d. Genome-wide association results for simulated phenotype %d\n\n'
            'Each panel is identified by its method name. The ten configurations are '
            'BOLT-LMM, REGENIE_LOCO, REGENIE_woLOCO, fastGWA, FarmCPU (NYC), '
            'FarmCPU (Kansas), BLINK (NYC), BLINK (Kansas), SLEMM and GCTA-MLMA. '
            'Every one of the 135,358 tested SNPs appears in each panel. Chromosomes '
            '1–18 share the same concatenated coordinates and alternate gray shades; '
            'a subset of chromosome numbers is labelled. Significant SNPs '
            '(P ≤ 5 × 10⁻⁸) are highlighted in the method color. The dashed gray line '
            'marks the significance threshold. All 100 panels in Figures S1–S10 '
            'share the ordinate range 0–110, matching Figure 4. Every observed '
            'value is displayed without truncation; method names appear above the axes. '
            'Significant SLEMM and GCTA-MLMA SNPs use larger circles for visibility. '
            'REGENIE LOG10P values are read directly; other values use −log₁₀(P).\n' % (ph,ph))


def main():
    global OUT, BENCH, GCTA
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark', type=Path, help='Authorized restricted benchmark root')
    parser.add_argument('--gcta-directory', type=Path, help='Completed restricted GCTA MLMA outputs')
    parser.add_argument('--source-cache', type=Path, help='Authorized local directory of 100 exported phenoXX_METHOD.tsv.gz tables')
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--counts', type=Path, default=DATA_ROOT / 'manhattan_counts.tsv')
    parser.add_argument('--main-only', action='store_true')
    parser.add_argument('--supplement-only', action='store_true')
    parser.add_argument('--phenotype', type=int, choices=range(1, 11), help='Render one supplementary phenotype only')
    args = parser.parse_args()
    if args.main_only and (args.supplement_only or args.phenotype):
        parser.error('--main-only cannot be combined with a supplementary-only selection')
    if args.source_cache is None and (args.benchmark is None or args.gcta_directory is None):
        parser.error('provide --source-cache or both --benchmark and --gcta-directory; variant-level inputs are restricted')
    OUT, BENCH, GCTA = args.output.resolve(), args.benchmark, args.gcta_directory
    for directory in ['audit', 'captions', 'supplement/captions', 'pdfs', 'png', 'supplement/pdfs', 'supplement/png']:
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    configure_style()
    expected = pd.read_csv(args.counts, sep='\t').set_index(['configuration', 'phenotype'])
    phenotypes = [args.phenotype] if args.phenotype else list(range(1, 11))
    data, records, inputs = {}, [], []
    for phenotype in phenotypes:
        for method in METHODS:
            if args.source_cache is None:
                frame, source, islog = load(method, phenotype)
            else:
                source = args.source_cache / ('pheno%02d_%s.tsv.gz' % (phenotype, method.replace('-', '_')))
                frame = pd.read_csv(source, sep='\t')
                islog = method.startswith('REGENIE')
                assert len(frame) == 135358 and not frame[['chr', 'position']].duplicated().any()
                assert np.isfinite(frame.native_value).all()
                assert np.array_equal(frame.significant.to_numpy(), significant(frame.native_value, islog))
                score = frame.native_value if islog else -np.log10(np.maximum(frame.native_value, 1e-300))
                assert np.allclose(frame.minus_log10p, score, rtol=1e-12, atol=1e-12)
                frame = frame.sort_values(['chr', 'position']).reset_index(drop=True)
            n = int(frame.significant.sum())
            detected = frame[frame.significant][['chr', 'position']].assign(window=lambda z: z.position // 500000)
            loci = len(detected.drop(columns='position').drop_duplicates())
            row = expected.loc[(method, phenotype)]
            assert n == int(row.significant_snps) and loci == int(row.loci), '{} {}'.format(method, phenotype)
            assert frame.minus_log10p.max() <= 110, 'Observed peak exceeds shared ordinate range'
            data[method, phenotype] = frame
            records.append(dict(configuration=method, phenotype=phenotype, n_tested=len(frame), significant_snps=n, loci=loci,
                                max_minus_log10p=float(frame.minus_log10p.max())))
            inputs.append(dict(configuration=method, phenotype=phenotype, filename=source.name, sha256=sha(source)))
        print('Validated phenotype {}'.format(phenotype), flush=True)
    reference = data[METHODS[0], phenotypes[0]]
    for frame in data.values():
        assert np.array_equal(frame[['chr', 'position']].values, reference[['chr', 'position']].values)
    lengths = reference.groupby('chr').position.max().sort_index()
    offsets = lengths.cumsum().shift(fill_value=0)
    centers = offsets + lengths / 2
    xmax = float((offsets.iloc[-1] + lengths.iloc[-1]) / 1e6)
    audit = []
    if not args.supplement_only and args.phenotype is None:
        main_figure(data, offsets, centers, xmax, audit)
    if not args.main_only:
        for phenotype in phenotypes:
            supplement(data, phenotype, offsets, centers, xmax, audit)
    write_captions(main_only=args.main_only, supplement_only=args.supplement_only or args.phenotype is not None)
    report = dict(status='PASS', significance='P <= 5e-8 on native reporting scale',
                  n_configurations=10, n_phenotypes=len(phenotypes),
                  all_counts_match_public_aggregates=True, shared_ordinate_range=[0, 110],
                  source_files=inputs, panel_counts=records, figures=audit,
                  source_sha256=sha(Path(__file__)))
    (OUT / 'audit/manhattan_validation.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
