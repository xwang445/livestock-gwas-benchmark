#!/usr/bin/env python3
"""Read-only physical stroke, text and raster audit of the 28 final PDFs.

Reads the actual Matplotlib PDF content streams, including Form XObjects.
Only painted strokes count as lines; fill-only geometry with width zero does
not. Graphics state, stroke opacity, affine transforms and q/Q are respected.
"""
import copy
import hashlib
import importlib.util
import json
import re
import subprocess
import zlib
from collections import Counter
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2] / 'figures/submission_2026-11-12'
SPEC = importlib.util.spec_from_file_location('submission_verifier', str(Path(__file__).resolve().with_name('verify_submission_figures.py')))
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


def tokens(data):
    i = 0
    length = len(data)
    separators = b'\x00\t\n\r\f ()<>[]{}/%'
    while i < length:
        c = data[i]
        if c in b'\x00\t\n\r\f ':
            i += 1
            continue
        if c == 37:
            stop = data.find(b'\n', i)
            i = length if stop < 0 else stop+1
            continue
        if c == 40:
            level = 1
            i += 1
            while i < length and level:
                if data[i] == 92:
                    i += 2
                    continue
                if data[i] == 40:
                    level += 1
                elif data[i] == 41:
                    level -= 1
                i += 1
            yield 'STRING'
            continue
        if c == 60 and data[i:i+2] != b'<<':
            stop = data.find(b'>', i+1)
            i = length if stop < 0 else stop+1
            yield 'STRING'
            continue
        if c in b'<>[]{}':
            yield chr(c)
            i += 1
            continue
        start = i
        if c == 47:
            i += 1
        while i < length and data[i] not in separators:
            i += 1
        if i == start:
            i += 1
        word = data[start:i].decode('latin1')
        try:
            yield float(word)
        except ValueError:
            yield word


def inspect_streams(pdf):
    raw = pdf.read_bytes()
    objects = {int(m.group(1)): m.group(2) for m in re.finditer(rb'(\d+) 0 obj\b(.*?)endobj', raw, re.S)}
    streams = {}
    for identifier, body in objects.items():
        hit = re.search(rb'\bstream\r?\n', body)
        if hit:
            # The final compressed byte may itself be CR; honor /Length rather
            # than trimming a possible CRLF before the endstream delimiter.
            length = re.search(rb'/Length\s+(\d+)(\s+0\s+R)?', body[:hit.start()])
            if not length:
                raise ValueError('Missing declared stream length')
            count = int(length.group(1))
            if length.group(2):
                count = int(objects[count].strip())
            value = body[hit.end():hit.end()+count]
            if len(value) != count:
                raise ValueError('Truncated PDF object body')
            if b'/FlateDecode' in body[:hit.start()]:
                value = zlib.decompress(value)
            streams[identifier] = value
    forms = {}
    for hit in re.finditer(rb'/([A-Za-z0-9_.-]+)\s+(\d+) 0 R', raw):
        name, identifier = hit.group(1).decode('latin1'), int(hit.group(2))
        if b'/Subtype /Form' in objects.get(identifier, b''):
            forms['/'+name] = identifier
    alphas = {}
    for hit in re.finditer(rb'/([A-Za-z0-9_.-]+)\s*<<([^<>]*?/Type /ExtGState[^<>]*?)>>', raw, re.S):
        body = hit.group(2)
        ca = re.search(rb'/CA\s+([\d.]+)', body)
        if ca:
            alphas['/'+hit.group(1).decode('latin1')] = float(ca.group(1))
    pages = [body for body in objects.values() if re.search(rb'/Type /Page\b', body)]
    content_ids = []
    for body in pages:
        hit = re.search(rb'/Contents\s+(\d+) 0 R', body)
        if not hit:
            raise ValueError('Unsupported page content representation')
        content_ids.append(int(hit.group(1)))
    sizes, declarations, stroke_widths, zero_strokes, ignored_alpha, form_calls = [], [], [], [], [], []
    initial = dict(width=1., alpha=1., matrix=np.eye(3), color=(0.,0.,0.))

    def walk(identifier, inherited, depth=0):
        if depth > 20:
            raise ValueError('Recursive Form XObject')
        state = copy.deepcopy(inherited)
        body = objects[identifier]
        transform = re.search(rb'/Matrix\s*\[([^]]+)\]', body)
        if transform:
            a,b,c,d,e,f = map(float,transform.group(1).split())
            state['matrix'] = state['matrix'].dot(np.array([[a,c,e],[b,d,f],[0,0,1.]]))
        stack, operands, path = [], [], False
        for token in tokens(streams[identifier]):
            if isinstance(token,float) or token.startswith('/'):
                operands.append(token)
                continue
            if token == 'q':
                stack.append(copy.deepcopy(state))
            elif token == 'Q':
                if not stack:
                    raise ValueError('Unbalanced PDF graphics state')
                state = stack.pop()
            elif token == 'w':
                state['width'] = float(operands[-1])
                declarations.append(state['width'])
            elif token == 'cm':
                a,b,c,d,e,f = operands[-6:]
                state['matrix'] = state['matrix'].dot(np.array([[a,c,e],[b,d,f],[0,0,1.]]))
            elif token == 'G':
                state['color'] = (operands[-1],)*3
            elif token == 'RG':
                state['color'] = tuple(operands[-3:])
            elif token == 'gs':
                name = operands[-1]
                if name not in alphas:
                    raise ValueError('Unknown stroke-opacity graphics state '+str(name))
                state['alpha'] = alphas[name]
            elif token == 'Tf':
                sizes.append(float(operands[-1]))
            elif token in ['m','l','c','v','y','h','re']:
                path = True
            elif token in ['S','s','B','B*','b','b*']:
                if path:
                    scale = float(np.linalg.svd(state['matrix'][:2,:2],compute_uv=False).min())
                    effective = state['width']*scale
                    detail = dict(object_id=identifier,operator=token,declared_width_pt=state['width'],
                        transform_minimum_scale=scale,effective_minimum_width_pt=effective,
                        stroke_alpha=state['alpha'],stroke_rgb=list(state['color']))
                    if state['alpha'] <= 0:
                        ignored_alpha.append(detail)
                    elif state['width'] == 0:
                        zero_strokes.append(detail)
                    else:
                        stroke_widths.append(detail)
                path = False
            elif token in ['f','F','f*','n']:
                path = False
            elif token == 'Do' and operands and operands[-1] in forms:
                child = forms[operands[-1]]
                form_calls.append(child)
                walk(child,state,depth+1)
            operands = []
        if stack:
            raise ValueError('Unbalanced PDF q/Q after content stream')

    for identifier in content_ids:
        walk(identifier,initial)
    return dict(actual_font_sizes_pt=sorted(set(sizes)),
        declared_linewidths_pt=sorted(set(declarations)),
        visible_stroke_count=len(stroke_widths),
        minimum_visible_stroke_pt=min([r['effective_minimum_width_pt'] for r in stroke_widths]) if stroke_widths else None,
        visible_stroke_width_histogram=dict(Counter(str(round(r['effective_minimum_width_pt'],8)) for r in stroke_widths)),
        visible_strokes_at_or_below_025pt=[r for r in stroke_widths if r['effective_minimum_width_pt'] <= .25],
        visible_zero_width_strokes=zero_strokes,
        transparent_stroke_count=len(ignored_alpha),
        zero_width_declarations=sum(v==0 for v in declarations),
        zero_width_note='Width-zero declarations for fill-only objects are not visible strokes; a width-zero painted stroke would fail.',
        form_xobject_calls=len(form_calls),parsed_page_content_streams=len(content_ids))


def main():
    previous=json.loads((ROOT/'audit/export_checks.json').read_text())
    old={r['figure']:r for r in previous['records']}
    records=[]
    for item in VERIFY.expected_inventory(ROOT):
        pdf=item['pdf']
        r=inspect_streams(pdf)
        font_lines=subprocess.check_output(['pdffonts',str(pdf)]).decode().splitlines()[2:]
        fonts=[line.split() for line in font_lines if line.strip()]
        image_lines=subprocess.check_output(['pdfimages','-list',str(pdf)]).decode().splitlines()[2:]
        rasters=[line.split() for line in image_lines if line.strip()]
        small=[v for v in r['actual_font_sizes_pt'] if v<8]
        fragments=old[item['figure']]['pdf']['small_font_fragments']
        math_allowed=(not small) or (bool(fragments) and all(f['math_exempt'] for f in fragments))
        with Image.open(str(item['png'])) as im:
            dpi=list(im.info.get('dpi',[]))
        checks=dict(visible_lines_strictly_above_025pt=not r['visible_strokes_at_or_below_025pt'],
            no_visible_zero_width_strokes=not r['visible_zero_width_strokes'],
            embedded_fonts=bool(fonts) and all(row[-5]=='yes' for row in fonts),
            no_type3_fonts=all('Type3' not in ''.join(row) for row in fonts),
            base_fonts_at_least_8pt=math_allowed and min(v for v in r['actual_font_sizes_pt'] if v>=8)>=8,
            raster_points_at_least_300dpi=all(float(row[12])>=300 and float(row[13])>=300 for row in rasters),
            png_300dpi=bool(dpi) and all(abs(v-300)<.2 for v in dpi))
        r.update(figure=item['figure'],pdf=str(pdf.relative_to(ROOT)),sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
            mathematical_script_sizes_below_8pt=small,mathematical_script_context_verified=math_allowed,
            minimum_base_font_pt=min(v for v in r['actual_font_sizes_pt'] if v>=8),
            raster_image_count=len(rasters),minimum_raster_dpi=min([float(row[12]) for row in rasters]) if rasters else None,
            png_dpi=dpi,checks=checks,status='PASS' if all(checks.values()) else 'FAIL')
        records.append(r)
        print(item['figure'],r['status'],'min_stroke',r['minimum_visible_stroke_pt'],'zero_strokes',len(r['visible_zero_width_strokes']),flush=True)
    output=dict(status='PASS' if all(r['status']=='PASS' for r in records) else 'FAIL',
        expected_figures=28,checked_figures=len(records),passed_figures=sum(r['status']=='PASS' for r in records),
        methodology='Actual PDF graphics-state/content-operator inspection, including Form XObjects, q/Q, stroke opacity and affine scale. Only painted strokes count; filled markers without edges are not strokes. Poppler confirms embedding/raster PPI; context-qualified math scripts are checked against the independent exported-artwork audit.',
        linewidth_limit_pt=.25,comparison='strictly greater than',records=records)
    (ROOT/'audit/all_saved_pdf_linewidths.json').write_text(json.dumps(output,indent=2)+'\n')
    return 0 if output['status']=='PASS' else 1


if __name__=='__main__':
    raise SystemExit(main())
