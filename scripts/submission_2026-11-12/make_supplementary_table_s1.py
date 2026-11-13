#!/usr/bin/env python3
"""Build editable and vector versions of the current Supplementary Table S1.

The primary table contains four columns of pooled regional percentages.
Counts and phenotype-level statistics are supplied in separate supporting
TSVs. The canonical analysis tables are read without edits.
"""
import argparse
import hashlib
import json
import math
import os
import zipfile
from datetime import datetime, timezone
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gwas_submission_matplotlib_20261112')
import numpy as np
import pandas as pd
from xml.etree import ElementTree as etree
from matplotlib.font_manager import FontProperties
from scipy.stats import t

from figure_style import ROOT, DATA_ROOT, METHOD_ORDER, LABELS, configure, plt

TITLE = ('Distribution of detected loci across polygenic-effect windows, null windows, '
         'and central-null regions.')
TABLE_ORDER = ['BOLT-LMM', 'REGENIE_LOCO', 'REGENIE_WO_LOCO', 'fastGWA',
               'SLEMM', 'GCTA_MLMA', 'FarmCPU_NYC', 'FarmCPU_Kansas', 'BLINK_NYC', 'BLINK_Kansas']
NOTES = [
    'Loci are fixed, chromosome-specific, non-overlapping 500-kb windows containing at least '
    'one SNP with P ≤ 5 × 10⁻⁸.',
    'Percentages pool locus detections across ten simulated phenotypes. The denominator is all '
    'detected locus–phenotype pairs for each method; a window detected in two phenotypes counts '
    'twice. These are not averages of phenotype-specific percentages.',
    'Central-null loci are a subset of null loci, not an additional mutually exclusive category.',
    'Values are rounded to two decimal places. Full pooled counts and phenotype-level results '
    'are provided in the supporting TSVs.',
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rounded(value, decimals=1):
    quantum = Decimal(1).scaleb(-decimals)
    return format(Decimal(str(float(value))).quantize(quantum, rounding=ROUND_HALF_UP),
                  '.{}f'.format(decimals))


def prepare(tables):
    by = pd.read_csv(tables / 'Table_S1_by_phenotype.tsv', sep='\t')
    by = by[by.configuration.isin(METHOD_ORDER)].copy()
    require(len(by) == 100 and not by.duplicated(['configuration', 'phenotype']).any(),
            'Expected one record for each of ten configurations and ten phenotypes.')
    by = by.set_index(['configuration', 'phenotype']).reindex(
        pd.MultiIndex.from_product([TABLE_ORDER, range(1, 11)], names=['configuration', 'phenotype'])
    ).reset_index()
    require(by.loci.notna().all(), 'Missing phenotype results.')
    require(by.central_null_windows.eq(1480).all(), 'Unexpected central-null denominator.')
    require(by.effect_loci.add(by.null_loci).eq(by.loci).all(), 'Effect/null counts do not partition loci.')
    require(by.central_null_loci.le(by.null_loci).all(), 'Central-null count exceeds null count.')
    require(np.allclose(by.central_null_spillover_rate, by.central_null_loci / 1480,
                        rtol=0, atol=1e-14), 'Spillover-rate mismatch.')

    background = pd.read_csv(tables / 'Table_S1_genomic_background.tsv', sep='\t', keep_default_na=False)
    summary = pd.read_csv(tables / 'Table_S1_locus_distribution.tsv', sep='\t').set_index('configuration')
    distributions, spillovers = [], []
    critical = float(t.ppf(.975, 9))
    for method in TABLE_ORDER:
        z = by[by.configuration.eq(method)]
        require(len(z) == 10, 'Every spillover mean must use all ten phenotypes.')
        total = int(z.loci.sum())
        effect = int(z.effect_loci.sum())
        null = int(z.null_loci.sum())
        central = int(z.central_null_loci.sum())
        require(total == int(summary.loc[method, 'pooled_loci']), 'Canonical pooled total mismatch.')
        record = dict(configuration=method, display_method=LABELS[method], n_phenotypes=10,
                      phenotypes_with_detected_loci=int(z.loci.gt(0).sum()),
                      pooled_loci=total, pooled_effect_loci=effect,
                      pooled_effect_percent=100 * effect / total,
                      pooled_null_loci=null, pooled_null_percent=100 * null / total,
                      pooled_central_null_loci=central, pooled_central_null_percent=100 * central / total)
        for field, canonical in [('pooled_effect_percent', 'pooled_effect_percent'),
                                 ('pooled_null_percent', 'pooled_null_percent'),
                                 ('pooled_central_null_percent', 'pooled_central_null_percent')]:
            require(np.isclose(record[field], summary.loc[method, canonical],
                               rtol=0, atol=1e-10), 'Canonical pooled fraction mismatch.')
        distributions.append(record)
        rates = z.central_null_loci.to_numpy(dtype=float) * 100 / 1480
        mean, sd = float(rates.mean()), float(rates.std(ddof=1))
        half_width = critical * sd / math.sqrt(10)
        spillovers.append(dict(configuration=method, display_method=LABELS[method], n_phenotypes=10,
                               central_null_windows_per_phenotype=1480,
                               central_null_window_phenotype_opportunities=14800,
                               central_null_loci_total=central,
                               central_null_loci_mean=float(z.central_null_loci.mean()),
                               central_null_loci_sd=float(z.central_null_loci.std(ddof=1)),
                               spillover_mean_percent=mean, spillover_sd_percent=sd,
                               spillover_ci95_lower_percent=mean - half_width,
                               spillover_ci95_upper_percent=mean + half_width,
                               t_df=9, t_critical=critical))
    distribution, spillover = pd.DataFrame(distributions), pd.DataFrame(spillovers)
    require(distribution.set_index('configuration').loc['BOLT-LMM', 'pooled_null_loci'] == 9592,
            'BOLT-LMM pooled null count must be 9,592.')
    require(distribution.set_index('configuration').loc['REGENIE_LOCO', 'pooled_null_loci'] == 9392,
            'REGENIE-LOCO pooled null count must be 9,392.')
    for method, total, effect in [('SLEMM', 16, 15), ('GCTA_MLMA', 15, 14)]:
        row = distribution.set_index('configuration').loc[method]
        require((row.pooled_loci, row.pooled_effect_loci, row.pooled_null_loci,
                 row.pooled_central_null_loci) == (total, effect, 1, 0), 'Current sparse-method inventory mismatch.')
    bg = background.set_index('region')
    require(bg.n_windows.to_dict() == {'effect': 1547, 'null': 2991, 'central_null': 1480},
            'Eligible-window inventory mismatch.')
    require(bg.tested_windows.to_dict() == {'effect': 1541, 'null': 2989, 'central_null': 1480},
            'Tested-window inventory mismatch.')
    require(background.eligible_windows.eq(4538).all() and background.total_tested_windows.eq(4530).all(),
            'Background total mismatch.')
    return by, distribution, spillover, background


def display_rows(distribution, spillover):
    rows = []
    for r in distribution.itertuples():
        rows.append([r.display_method, rounded(r.pooled_effect_percent, 2),
                     rounded(r.pooled_null_percent, 2), rounded(r.pooled_central_null_percent, 2)])
    return rows


W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
etree.register_namespace('w', W)
etree.register_namespace('r', R)


def element(name, attrs=None, text=None):
    node = etree.Element('{%s}%s' % (W, name))
    for key, value in (attrs or {}).items():
        node.set('{%s}%s' % (W, key), str(value))
    if text is not None:
        node.text = text
    return node


def paragraph(text, size=20, bold=False, after=70, before=0, color='222222', keep=False, align='left'):
    node = element('p')
    pp = element('pPr')
    pp.append(element('spacing', {'before': before, 'after': after, 'line': 250, 'lineRule': 'auto'}))
    pp.append(element('jc', {'val': align}))
    if keep:
        pp.append(element('keepNext'))
    node.append(pp)
    for index, line in enumerate(text.split('\n')):
        run = element('r')
        rp = element('rPr')
        rp.append(element('rFonts', {'ascii': 'Arial', 'hAnsi': 'Arial', 'eastAsia': 'Arial', 'cs': 'Arial'}))
        rp.append(element('sz', {'val': size}))
        rp.append(element('color', {'val': color}))
        if bold:
            rp.append(element('b'))
        run.append(rp)
        if index:
            run.append(element('br'))
        tnode = element('t', text=line)
        tnode.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
        run.append(tnode)
        node.append(run)
    return node


def word_table(headers, rows, widths_mm):
    table = element('tbl')
    pr = element('tblPr')
    pr.append(element('tblW', {'w': round(sum(widths_mm) / 25.4 * 1440), 'type': 'dxa'}))
    pr.append(element('tblLayout', {'type': 'fixed'}))
    borders = element('tblBorders')
    for edge in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
        borders.append(element(edge, {'val': 'nil'}))
    pr.append(borders)
    margins = element('tblCellMar')
    for side, value in [('top', 36), ('bottom', 36), ('left', 75), ('right', 75)]:
        margins.append(element(side, {'w': value, 'type': 'dxa'}))
    pr.append(margins)
    table.append(pr)
    grid = element('tblGrid')
    for width in widths_mm:
        grid.append(element('gridCol', {'w': round(width / 25.4 * 1440)}))
    table.append(grid)
    for row_index, values in enumerate([headers] + rows):
        tr = element('tr')
        trpr = element('trPr')
        trpr.append(element('cantSplit'))
        if row_index == 0:
            trpr.append(element('tblHeader'))
        tr.append(trpr)
        for col_index, (value, width) in enumerate(zip(values, widths_mm)):
            cell = element('tc')
            cp = element('tcPr')
            cp.append(element('tcW', {'w': round(width / 25.4 * 1440), 'type': 'dxa'}))
            cp.append(element('vAlign', {'val': 'center'}))
            if row_index == 0 or row_index % 2 == 0:
                cp.append(element('shd', {'fill': 'EDEFF2' if row_index == 0 else 'F7F8FA', 'val': 'clear'}))
            cb = element('tcBorders')
            if row_index == 0:
                for side in ['top', 'bottom']:
                    cb.append(element(side, {'val': 'single', 'sz': 5, 'color': '555555'}))
            elif row_index == len(rows):
                cb.append(element('bottom', {'val': 'single', 'sz': 5, 'color': '555555'}))
            cp.append(cb)
            cell.append(cp)
            cell.append(paragraph(str(value), size=18 if row_index == 0 else 19,
                                  bold=row_index == 0, after=0, keep=row_index == 0,
                                  align='left' if col_index == 0 else 'center'))
            tr.append(cell)
        table.append(tr)
    return table


def xml_bytes(node):
    return etree.tostring(node, encoding='UTF-8', xml_declaration=True)


def write_docx(path, headers, rows):
    document = etree.Element('{%s}document' % W)
    body = element('body')
    document.append(body)
    body.append(paragraph('Supplementary Table S1. ' + TITLE, size=24, bold=True,
                          color='19354A', after=100, keep=True))
    body.append(word_table(headers, rows, [50, 40, 40, 40]))
    body.append(paragraph('Notes', size=19, bold=True, before=130, after=45, keep=True))
    for note in NOTES:
        body.append(paragraph(note, size=18, after=55))
    section = element('sectPr')
    section.append(element('pgSz', {'w': 11906, 'h': 16838}))
    section.append(element('pgMar', {'top': 850, 'bottom': 850, 'left': 1134, 'right': 1134,
                                    'header': 360, 'footer': 360, 'gutter': 0}))
    body.append(section)
    styles = etree.Element('{%s}styles' % W)
    default = element('style', {'type': 'paragraph', 'default': '1', 'styleId': 'Normal'})
    default.append(element('name', {'val': 'Normal'}))
    styles.append(default)
    types = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
</Types>'''
    rels = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>'''
    doc_rels = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>'''
    with zipfile.ZipFile(str(path), 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('[Content_Types].xml', types)
        archive.writestr('_rels/.rels', rels)
        archive.writestr('word/document.xml', xml_bytes(document))
        archive.writestr('word/styles.xml', xml_bytes(styles))
        archive.writestr('word/_rels/document.xml.rels', doc_rels)
    with zipfile.ZipFile(str(path)) as archive:
        require(archive.testzip() is None, 'DOCX ZIP checksum failed.')
        root = etree.fromstring(archive.read('word/document.xml'))
        require(len(root.findall('.//{%s}tbl' % W)) == 1, 'DOCX must contain one editable table.')
        require(len(root.findall('.//{%s}tr' % W)) == 11, 'DOCX row inventory mismatch.')
        require(all(len(row.findall('{%s}tc' % W)) == 4 for row in root.findall('.//{%s}tr' % W)),
                'DOCX must have exactly four columns in every row.')


def pdf_table(ax, top, headers, rows, widths, left=20):
    x = left
    header_height, row_height = 10, 5.7
    total = sum(widths)
    ax.add_patch(plt.Rectangle((left, top - header_height), total, header_height,
                               facecolor='#EDEFF2', edgecolor='none', linewidth=.35))
    for index, (header, width) in enumerate(zip(headers, widths)):
        ax.text(x + (1.7 if index == 0 else width / 2), top - header_height / 2,
                header, ha='left' if index == 0 else 'center', va='center',
                fontsize=9, fontweight='bold', linespacing=1.18)
        x += width
    ax.plot([left, left + total], [top, top], color='#555555', linewidth=.65)
    ax.plot([left, left + total], [top - header_height, top - header_height],
            color='#555555', linewidth=.5)
    y = top - header_height
    for row_index, row in enumerate(rows):
        if row_index % 2:
            ax.add_patch(plt.Rectangle((left, y - row_height), total, row_height,
                                       facecolor='#F7F8FA', edgecolor='none', linewidth=.35))
        x = left
        for index, (value, width) in enumerate(zip(row, widths)):
            ax.text(x + (1.7 if index == 0 else width / 2), y - row_height / 2,
                    value, ha='left' if index == 0 else 'center', va='center', fontsize=10)
            x += width
        y -= row_height
    ax.plot([left, left + total], [y, y], color='#555555', linewidth=.65)
    return y


def write_pdf(path, png, headers, rows):
    configure(10)
    fig = plt.figure(figsize=(210 / 25.4, 190 / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 210), ylim=(0, 190))
    ax.axis('off')
    ax.text(20, 166, 'Supplementary Table S1', fontsize=12, fontweight='bold', color='#19354A')
    for index, line in enumerate(['Distribution of detected loci across polygenic-effect windows,',
                                  'null windows, and central-null regions.']):
        ax.text(20, 158.5 - 5.5 * index, line, fontsize=11, fontweight='bold')
    bottom = pdf_table(ax, 144, headers, rows, [50, 40, 40, 40])
    ax.text(20, bottom - 9, 'Notes', fontsize=9.2, fontweight='bold')
    y = bottom - 14
    renderer = fig.canvas.get_renderer()
    note_font = FontProperties(family='DejaVu Sans', size=9)
    maximum_width_px = 169.5 / 25.4 * fig.dpi
    for note in NOTES:
        lines, current = [], ''
        for word in note.split():
            candidate = current + (' ' if current else '') + word
            width = renderer.get_text_width_height_descent(candidate, note_font, False)[0]
            if current and width > maximum_width_px:
                lines.append(current)
                current = word
            else:
                current = candidate
        if current:
            lines.append(current)
        for line in lines:
            ax.text(20, y, line, fontsize=9, ha='left', va='top')
            y -= 3.85
        y -= 1.25
    require(y > 20, 'PDF notes exceed the bottom content margin.')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for text in fig.findobj(plt.Text):
        if not text.get_visible() or not text.get_text():
            continue
        b = text.get_window_extent(renderer)
        if b.x0 < 19.8 / 25.4 * fig.dpi or b.y0 < 0 or \
                b.x1 > 190.2 / 25.4 * fig.dpi or b.y1 > fig.bbox.height:
            outside.append(text.get_text())
    require(not outside, 'Text exceeds PDF page: ' + repr(outside))
    fig.savefig(str(path), format='pdf', metadata={'Title': 'Supplementary Table S1',
                                                  'Creator': 'Canonical benchmark tables; vector Matplotlib export'})
    fig.savefig(str(png), dpi=300)
    plt.close(fig)
    return dict(page_width_mm=210, page_height_mm=190, minimum_font_pt=9,
                text_within_page=True, text_within_170_mm_content_column=True,
                pdf_fonttype=42, numerical_table_text_is_vector=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path,
                        default=DATA_ROOT / 'tables')
    parser.add_argument('--output', type=Path, default=ROOT / 'supplement/tables')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    audit = output / 'audit'
    audit.mkdir(parents=True, exist_ok=True)
    source_names = ['Table_S1_by_phenotype.tsv', 'Table_S1_locus_distribution.tsv',
                    'Table_S1_genomic_background.tsv']
    source_hashes = {name: sha(args.tables / name) for name in source_names}
    by, distribution, spillover, background = prepare(args.tables)
    merged = distribution.merge(spillover.drop(columns=['display_method', 'n_phenotypes']), on='configuration',
                                validate='one_to_one')
    merged.to_csv(output / 'Table_S1_supporting_statistics.tsv', sep='\t', index=False, float_format='%.15g')
    distribution.to_csv(output / 'Table_S1_locus_distribution.tsv', sep='\t', index=False, float_format='%.15g')
    spillover.to_csv(output / 'Table_S1_central_null_spillover.tsv', sep='\t', index=False, float_format='%.15g')
    fields = ['configuration', 'phenotype', 'loci', 'effect_loci', 'null_loci', 'central_null_loci',
              'central_null_windows', 'central_null_spillover_rate']
    detail = by[fields].copy()
    detail.insert(1, 'display_method', detail.configuration.map(LABELS))
    detail['effect_percent_of_detected'] = 100 * detail.effect_loci / detail.loci.replace(0, np.nan)
    detail['null_percent_of_detected'] = 100 * detail.null_loci / detail.loci.replace(0, np.nan)
    detail['central_null_percent_of_detected'] = 100 * detail.central_null_loci / detail.loci.replace(0, np.nan)
    detail['central_null_spillover_percent'] = detail.central_null_loci * 100 / detail.central_null_windows
    detail.to_csv(output / 'Table_S1_by_phenotype.tsv', sep='\t', index=False,
                  float_format='%.15g', na_rep='NA')
    background.to_csv(output / 'Table_S1_genomic_background.tsv', sep='\t', index=False, float_format='%.15g')
    headers = ['Method', 'Effect %', 'Null %', 'Central-null %']
    rows = display_rows(distribution, spillover)
    pd.DataFrame(rows, columns=headers).to_csv(output / 'Table_S1.tsv', sep='\t', index=False)
    write_docx(output / 'Table_S1.docx', headers, rows)
    layout = write_pdf(output / 'Table_S1.pdf', output / 'Table_S1_preview.png', headers, rows)
    caption_text = ('Supplementary Table S1. ' + TITLE + '\n\n' +
                    '\n\n'.join(NOTES) + '\n')
    (output / 'Table_S1_caption.txt').write_text(caption_text)
    pd.DataFrame(rows, columns=headers).to_csv(output / 'Table_S1_display_locus_distribution.tsv',
                                              sep='\t', index=False)
    # The primary table now has distribution percentages only. Retain the
    # full-precision spillover TSV as supporting analysis, but remove its
    # obsolete display table so no old primary layout can be mistaken for S1.
    obsolete = output / 'Table_S1_display_spillover.tsv'
    if obsolete.exists():
        obsolete.unlink()
    require(source_hashes == {name: sha(args.tables / name) for name in source_names},
            'A canonical source changed during generation.')
    source_manifest = pd.DataFrame([dict(path=name,
                                         sha256=value, bytes=(args.tables / name).stat().st_size)
                                    for name, value in source_hashes.items()])
    source_manifest.to_csv(audit / 'table_s1_source_manifest.tsv', sep='\t', index=False)
    validation = dict(status='PASS', table='Supplementary Table S1',
                      created_utc=datetime.now(timezone.utc).isoformat(),
                      n_configurations=10, n_phenotypes=10, n_phenotype_rows=len(by),
                      all_pooled_counts_recomputed_from_phenotype_aggregates=True,
                      all_pooled_fractions_match_canonical_summary=True,
                      all_spillover_means_recomputed_from_ten_phenotypes=True,
                      primary_columns=headers, primary_table_shape_including_header=[11, 4],
                      primary_rows=rows, primary_row_order=TABLE_ORDER,
                      primary_table_contains_spillover_section=False,
                      supporting_spillover_CI_definition='mean ± t(0.975,9) × sample SD / sqrt(10), percentage points',
                      supporting_t_critical=float(t.ppf(.975, 9)),
                      all_ten_phenotypes_included_in_supporting_CIs=True,
                      locus_distribution_aggregation='pooled detections across ten phenotypes',
                      central_null_is_subset_of_null=True,
                      eligible_windows=4538, tested_SNP_containing_windows=4530,
                      central_null_windows_per_phenotype=1480,
                      current_SLEMM_counts=dict(total=16, effect=15, null=1, central_null=0),
                      current_GCTA_MLMA_counts=dict(total=15, effect=14, null=1, central_null=0),
                      pooled_statistics=distribution.to_dict(orient='records'),
                      display_rounding='two decimals, decimal ROUND_HALF_UP',
                      source_categories_preserved_with_keep_default_na_false=True,
                      source_hashes_unchanged_during_generation=True,
                      docx_contains_one_editable_table=True, docx_table_rows_including_header=11,
                      docx_table_columns=4, sections_have_no_panel_letters=True, pdf_layout=layout,
                      sources=source_manifest.to_dict(orient='records'), script_sha256=sha(__file__),
                      outputs=[dict(path=p.name, bytes=p.stat().st_size, sha256=sha(p))
                               for p in sorted(output.iterdir()) if p.is_file()])
    (audit / 'table_s1_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(dict(status='PASS', output=str(output),
                          audit=str(audit / 'table_s1_validation.json'),
                          pooled_loci=distribution[['configuration', 'pooled_loci']].to_dict(orient='records'),
                          primary_display=rows), indent=2))


if __name__ == '__main__':
    main()
