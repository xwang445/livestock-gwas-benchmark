#!/usr/bin/env python3
"""Independent exported-artwork checks for eight main and twenty supplements.

Requires Poppler and Pillow. Missing files are reported explicitly, rather than
being silently skipped. Titles and captions are checked outside the artwork.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

from PIL import Image, ImageChops, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2] / 'figures/submission_2026-11-12'
SIZE_LIMIT = 10000000
NUMBERED_FIGURE = re.compile(r'\b(?:Supplementary\s+)?Figure\s+S?\d+\b', re.I)
CAPTION_PREFIX = re.compile(r'^\s*(?:Supplementary\s+)?Figure\s+S?\d+[.:]?\s*', re.I)
ALIAS_PATTERNS = [r'\bGCTA MLMA\b', r'\bREGENIE-LOCO\b', r'\bREGENIE_WO_LOCO\b',
                  r'\bSLEMM_previous\b']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def command(argv):
    result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError('{}: {}'.format(' '.join(argv), result.stderr.decode('utf-8', errors='replace')))
    return result.stdout.decode('utf-8')


def pick(candidates):
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


def expected_inventory(root):
    records = []
    for number in range(1, 9):
        stem = 'Figure_{:02d}'.format(number)
        records.append(dict(
            figure=stem, category='main',
            pdf=root/'pdfs'/(stem+'.pdf'), png=root/'png'/(stem+'.png'),
            caption=pick([root/'captions'/(stem+'_caption.txt'),
                          root/'captions_draft'/(stem+'_caption.txt')]),
        ))
    for number in range(1, 21):
        stem = 'Figure_S{:02d}'.format(number)
        records.append(dict(
            figure=stem, category='supplement',
            pdf=pick([root/'supplement/pdfs'/(stem+'.pdf'), root/'supplement/pdf'/(stem+'.pdf')]),
            png=root/'supplement/png'/(stem+'.png'),
            caption=pick([root/'supplement/captions'/(stem+'_caption.txt'),
                          root/'captions_draft'/(stem+'_caption.txt')]),
        ))
    return records


def caption_check(path):
    text = path.read_text(encoding='utf-8').strip()
    unnumbered = CAPTION_PREFIX.sub('', text, count=1)
    parts = re.split(r'\.(?:\s+|$)|\n', unnumbered, maxsplit=1)
    title = parts[0].strip()
    legend = parts[1].strip() if len(parts) > 1 else ''
    return dict(title=title, title_words=len(title.split()), legend_words=len(legend.split()),
                title_pass=len(title.split()) <= 15,
                legend_pass=bool(legend) and len(legend.split()) <= 300,
                caption_total_words=len(text.split()),
                utf8_valid=True)


def math_fragment(text):
    """Recognize short superscript/subscript text; Poppler rounds 8.5 pt to 8."""
    value = ''.join(text.split())
    if not value:
        return True
    if re.fullmatch(r'[0-9⁰¹²³⁴⁵⁶⁷⁸⁹₀₁₂₃₄₅₆₇₈₉.,−+×/()]+', value) and len(value) <= 8:
        return True
    return value in {'GC', 'NCP', 'n', 'm', 'p', 'r', 'h', 'e', 'g', 'x', 'y', 't',
                     'α', 'β', 'χ', 'λ', 'σ', 'μ', 'i', 'j', 'k'}


def has_math_context(nodes, index, specs):
    """A small numeral in a heatmap is not automatically a mathematical script."""
    node = nodes[index]
    x, y = float(node.get('left')), float(node.get('top'))
    w, h = float(node.get('width')), float(node.get('height'))
    neighbors = nodes[max(0,index-2):index] + nodes[index+1:index+3]
    for neighbor in neighbors:
        if specs.get(neighbor.get('font'), 0) < 8:
            continue
        text = ''.join(neighbor.itertext()).strip()
        if not re.search(r'[χλ]|×10|log|\br\b|\bh\b|\bp\b|\bt\b|NCP', text):
            continue
        nx, ny = float(neighbor.get('left')), float(neighbor.get('top'))
        nw, nh = float(neighbor.get('width')), float(neighbor.get('height'))
        gap_x = max(x-(nx+nw), nx-(x+w), 0)
        gap_y = max(y-(ny+nh), ny-(y+h), 0)
        if gap_x <= 8 and gap_y <= 8:
            return True
    return False


def pdf_check(path):
    info = command(['pdfinfo', str(path)])
    fields = {}
    for line in info.splitlines():
        if ':' in line:
            key, value = line.split(':', 1)
            fields[key.strip()] = value.strip()
    dimensions = re.search(r'([0-9.]+)\s+x\s+([0-9.]+)\s+pts', fields.get('Page size', ''))
    if not dimensions:
        raise ValueError('No PDF page dimensions: ' + str(path))
    width_pt, height_pt = map(float, dimensions.groups())
    width_mm, height_mm = width_pt*25.4/72, height_pt*25.4/72
    fonts = command(['pdffonts', str(path)])
    font_lines = [line for line in fonts.splitlines()[2:] if line.strip()]
    embedding = [line.split()[-5] == 'yes' for line in font_lines]
    xml = command(['pdftohtml', '-xml', '-i', '-hidden', '-zoom', '1', '-stdout', str(path)])
    tree = ET.fromstring(xml)
    specs = {node.get('id'): float(node.get('size')) for node in tree.findall('.//fontspec')}
    small_fragments, base_sizes, all_sizes, outside = [], [], [], []
    for page in tree.findall('page'):
        page_width, page_height = float(page.get('width')), float(page.get('height'))
        text_nodes = page.findall('text')
        for index, node in enumerate(text_nodes):
            text = ''.join(node.itertext()).strip()
            if not text:
                continue
            size = specs.get(node.get('font'))
            if size is not None:
                all_sizes.append(size)
                math_exempt = size < 8 and math_fragment(text) and has_math_context(text_nodes, index, specs)
                if size < 8:
                    small_fragments.append(dict(size_pt=size, text=text, math_exempt=math_exempt))
                if not math_exempt:
                    base_sizes.append(size)
            left, top = float(node.get('left')), float(node.get('top'))
            width, height = float(node.get('width')), float(node.get('height'))
            # Integer XML coordinates can extend by about a point due to rounding.
            if left < -1 or top < -1 or left+width > page_width+1 or top+height > page_height+1:
                outside.append(dict(text=text, box=[left, top, width, height]))
    text = command(['pdftotext', '-layout', str(path), '-'])
    invalid_chars = sorted(set(char for char in text if char in '\ufffd\u25a1\x00'))
    bad_aliases = [pattern for pattern in ALIAS_PATTERNS if re.search(pattern, text)]
    record = dict(
        bytes=path.stat().st_size, sha256=sha256(path), pages=int(fields.get('Pages', 0)),
        width_mm=width_mm, height_mm=height_mm, rotation=int(fields.get('Page rot', 0)),
        embedded_fonts=bool(font_lines) and all(embedding), font_count=len(font_lines),
        fontspec_sizes_pt=sorted(set(all_sizes)),
        minimum_base_font_pt=min(base_sizes) if base_sizes else None,
        small_font_fragments=small_fragments, text_outside_page=outside,
        numbered_figure_titles_in_artwork=NUMBERED_FIGURE.findall(text),
        invalid_unicode_characters=invalid_chars, inconsistent_method_aliases=bad_aliases,
        text_extracted_characters=len(text),
    )
    record['checks'] = dict(
        file_size=record['bytes'] < SIZE_LIMIT,
        width=abs(width_mm-170) <= .2,
        height=0 < height_mm <= 190.01,
        single_page=record['pages'] == 1,
        correct_orientation=record['rotation'] == 0,
        embedded_fonts=record['embedded_fonts'],
        no_numbered_title_or_caption=not record['numbered_figure_titles_in_artwork'],
        unicode=not invalid_chars,
        method_names=not bad_aliases,
        base_font_size=record['minimum_base_font_pt'] is not None and record['minimum_base_font_pt'] >= 8,
        no_text_outside_page=not outside,
    )
    return record, info, fonts


def png_check(path, pdf):
    with Image.open(path) as image:
        width, height = image.size
        dpi = image.info.get('dpi')
        rgb = image.convert('RGB')
        bbox = ImageChops.difference(rgb, Image.new('RGB', rgb.size, 'white')).getbbox()
    expected_width = pdf['width_mm'] / 25.4 * 300
    expected_height = pdf['height_mm'] / 25.4 * 300
    dpi_ok = bool(dpi) and all(abs(float(v)-300) <= .2 for v in dpi)
    margin_fractions = None
    if bbox:
        margin_fractions = [bbox[0]/width, bbox[1]/height,
                            (width-bbox[2])/width, (height-bbox[3])/height]
    return dict(bytes=path.stat().st_size, sha256=sha256(path), width_px=width, height_px=height,
                dpi=dpi, content_bbox_px=bbox, outer_margin_fractions=margin_fractions,
                checks=dict(file_size=path.stat().st_size < SIZE_LIMIT,
                            dpi=dpi_ok,
                            dimensions_match_pdf=abs(width-expected_width) <= 2 and abs(height-expected_height) <= 2,
                            visible_artwork=bool(bbox),
                            close_outer_crop=bool(margin_fractions) and max(margin_fractions) <= .15))


def contact_sheet(records, destination, main_only=False):
    available = [record for record in records if Path(record['png_path']).exists()
                 and (not main_only or record['category'] == 'main')]
    if not available:
        return
    columns = 2 if main_only else 4
    box_width, box_height = (610, 685) if main_only else (360, 455)
    rows = (len(available)+columns-1)//columns
    sheet = Image.new('RGB', (columns*box_width, rows*box_height), 'white')
    draw = ImageDraw.Draw(sheet)
    font_path = '/usr/share/fonts/dejavu/DejaVuSans.ttf'
    try:
        font = ImageFont.truetype(font_path, 18 if main_only else 14)
    except OSError:
        font = ImageFont.load_default()
    for index, record in enumerate(available):
        x, y = (index % columns)*box_width, (index // columns)*box_height
        with Image.open(record['png_path']) as image:
            thumb = image.convert('RGB')
            thumb.thumbnail((box_width-20, box_height-45), Image.LANCZOS)
        draw.text((x+10, y+8), record['figure'], fill='#222222', font=font)
        sheet.paste(thumb, (x+(box_width-thumb.width)//2, y+35))
    sheet.save(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--allow-missing', action='store_true', help='Report incomplete state without a nonzero exit.')
    parser.add_argument('--contact-sheets', action='store_true')
    args = parser.parse_args()
    audit = args.root/'audit'
    audit.mkdir(parents=True, exist_ok=True)
    records, missing, failures = [], [], []
    for item in expected_inventory(args.root):
        stem = item['figure']
        record = dict(figure=stem, category=item['category'],
                      pdf_path=str(item['pdf'].relative_to(args.root)), png_path=str(item['png'].relative_to(args.root)),
                      caption_path=str(item['caption'].relative_to(args.root)), status='pending')
        absent = [str(item[k].relative_to(args.root)) for k in ['pdf', 'png', 'caption'] if not item[k].exists()]
        if absent:
            record['status'], record['missing'] = 'missing', absent
            missing.extend(absent)
            records.append(record)
            continue
        try:
            pdf, info, fonts = pdf_check(item['pdf'])
            png = png_check(item['png'], pdf)
            cap = caption_check(item['caption'])
            (audit/(stem+'_pdfinfo.txt')).write_text(info)
            (audit/(stem+'_pdffonts.txt')).write_text(fonts)
            record.update(pdf=pdf, png=png, caption=cap)
            failed = ['pdf.'+k for k,v in pdf['checks'].items() if not v]
            failed += ['png.'+k for k,v in png['checks'].items() if not v]
            failed += ['caption.'+k for k in ['title_pass', 'legend_pass'] if not cap[k]]
            record['status'] = 'fail' if failed else 'pass'
            record['failed_checks'] = failed
            failures.extend(stem+': '+value for value in failed)
        except Exception as error:
            record['status'], record['error'] = 'error', str(error)
            failures.append(stem+': '+str(error))
        records.append(record)
    complete = not missing and not failures
    payload = dict(status='PASS' if complete else ('INCOMPLETE' if missing and not failures else 'FAIL'),
                   expected_figures=28, main_figures=8, supplementary_figures=20,
                   checked_figures=sum(r['status'] in ['pass','fail'] for r in records),
                   passed_figures=sum(r['status']=='pass' for r in records),
                   missing=missing, failures=failures, records=records,
                   limits=dict(width_mm=170, width_tolerance_mm=.2, graphic_height_mm=190,
                               file_bytes=SIZE_LIMIT, png_dpi=300,
                               base_font_poppler_rounded_pt=8, title_words=15, legend_words=300),
                   font_note='Poppler XML at zoom 1 reports integer sizes; 8.5 pt base type is reported as 8. Short mathematical super/subscripts are exempt from the base-size check.',
                   source_sha256=sha256(__file__))
    (audit/'export_checks.json').write_text(json.dumps(payload, indent=2, ensure_ascii=False)+'\n')
    columns = ['figure','category','status','width_mm','height_mm','pdf_bytes','png_bytes','font_count',
               'minimum_base_font_pt','title_words','legend_words','failed_checks']
    with (audit/'export_checks.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for record in records:
            pdf, png, cap = record.get('pdf',{}), record.get('png',{}), record.get('caption',{})
            writer.writerow(dict(figure=record['figure'], category=record['category'], status=record['status'],
                                 width_mm=pdf.get('width_mm'), height_mm=pdf.get('height_mm'),
                                 pdf_bytes=pdf.get('bytes'), png_bytes=png.get('bytes'),
                                 font_count=pdf.get('font_count'), minimum_base_font_pt=pdf.get('minimum_base_font_pt'),
                                 title_words=cap.get('title_words'), legend_words=cap.get('legend_words'),
                                 failed_checks='; '.join(record.get('failed_checks',[]))))
    if args.contact_sheets:
        display_records = [dict(record, png_path=str(args.root / record['png_path'])) for record in records]
        contact_sheet(display_records, audit/'main_figures_contact.png', main_only=True)
        contact_sheet(display_records, audit/'all_figures_contact.png')
    print('{}: {} / 28 figures pass; {} missing files; {} failed checks.'.format(
        payload['status'], payload['passed_figures'], len(missing), len(failures)))
    for value in failures:
        print('FAIL:', value)
    for value in missing:
        print('MISSING:', value)
    return 0 if complete or (args.allow_missing and not failures) else 1


if __name__ == '__main__':
    sys.exit(main())
