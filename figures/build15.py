# -*- coding: utf-8 -*-
"""Build draft 15 (main manuscript) for Computational Materials Science.

Style carriers and equations come from draft_12.docx; tables 1 and 3-9 are taken from draft_14.docx and restyled
(no vertical rules, no shading); Table 2 is rebuilt so that its citations renumber with the text.
"""
import sys, copy, re, json
sys.path.insert(0, '/tmp/w/b15'); sys.path.insert(1, '/tmp/w/b14')
import docx
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph
from docx.table import Table
from lib import build_table, settext, set_cell
import t15_front as FRONT, t15_a, t15_b, t15_d
from refs15 import REFS
import edits15
_miss = edits15.apply((t15_a, t15_b, t15_d))
if _miss: raise SystemExit('EDIT NOT APPLIED: %s' % _miss)

SRC12 = '/mnt/user-data/outputs/draft_12.docx'
SRC14 = '/mnt/user-data/outputs/draft_14.docx'
OUT = '/tmp/w/b15/draft_15.docx'
XMLSP = '{http://www.w3.org/XML/1998/namespace}space'

doc = docx.Document(SRC12)
body = doc.element.body
P = list(doc.paragraphs)

def grab(i):
    el = copy.deepcopy(P[i]._p)
    for r in el.findall(qn('w:r')):
        if r.find(qn('w:drawing')) is not None:
            el.remove(r)
    return el

CAP_IDX = next(i for i, p in enumerate(P)
               if ''.join(r.text for r in p.runs).strip().startswith('Table 1. Composition matrix'))
CAR = {'title': grab(0), 'auth': grab(1), 'H1': grab(2), 'body': grab(3), 'H2': grab(18), 'cap': grab(CAP_IDX)}
EQ = {'(1)': copy.deepcopy(P[44]._p), '(2)': copy.deepcopy(P[51]._p),
      '(3)': copy.deepcopy(P[52]._p), '(4)': copy.deepcopy(P[57]._p)}
T14 = [copy.deepcopy(t._tbl) for t in docx.Document(SRC14).tables]

for child in list(body):
    if child.tag != qn('w:sectPr'):
        body.remove(child)

def _append(el):
    body.append(el)

def new_par(kind):
    el = copy.deepcopy(CAR[kind])
    for r in list(el):
        if r.tag in (qn('w:r'), qn('w:proofErr'), qn('w:bookmarkStart'), qn('w:bookmarkEnd'), qn('w:hyperlink')):
            el.remove(r)
    _append(el)
    return Paragraph(el, doc._body)

def put(kind, text, align=None):
    el = copy.deepcopy(CAR[kind]); _append(el)
    p = Paragraph(el, doc._body); settext(p, text)
    if align is not None:
        p.alignment = align
    return p

def put_runs(kind, runs, align=None, size=None):
    """runs: list of (text, flag) with flag in {False, True='sup', 'b'}"""
    p = put(kind, runs[0][0], align)
    base = p.runs[0]._r
    for r in list(p._p.findall(qn('w:r'))):
        p._p.remove(r)
    for text, flag in runs:
        r = copy.deepcopy(base)
        for t in r.findall(qn('w:t')):
            r.remove(t)
        rpr = r.find(qn('w:rPr'))
        if rpr is None:
            rpr = OxmlElement('w:rPr'); r.insert(0, rpr)
        if flag is True:
            va = OxmlElement('w:vertAlign'); va.set(qn('w:val'), 'superscript'); rpr.append(va)
        if flag == 'b':
            rpr.insert(0, OxmlElement('w:b'))
        if size:
            for tag in ('w:sz', 'w:szCs'):
                e = rpr.find(qn(tag))
                if e is None:
                    e = OxmlElement(tag); rpr.append(e)
                e.set(qn('w:val'), str(int(size * 2)))
        t = OxmlElement('w:t'); t.text = text; t.set(XMLSP, 'preserve'); r.append(t)
        p._p.append(r)
    return p

def keep_next(p):
    pPr = p._p.get_or_add_pPr()
    if pPr.find(qn('w:keepNext')) is None:
        pPr.append(OxmlElement('w:keepNext'))

def put_caption(spec):
    label, text = spec.split('|', 1)
    p = put_runs('cap', [(label + ' ', 'b'), (text, False)], AL.JUSTIFY)
    return p

def put_fig(path, width_in, capt):
    p = new_par('body'); p.alignment = AL.CENTER
    p.add_run().add_picture(path, width=Inches(width_in or 6.2))
    keep_next(p)
    put_caption(capt)

# ---------------------------------------------------------------- Table 2 (rebuilt, keyed citations)
MINUS = '−'
t2_rows = [
    ['Mg lattice constant, a', '3.209 Å', '3.21 Å (experiment) [@kittel]'],
    ['c/a ratio', '1.619', '1.623 (experiment) [@kittel]'],
    ['Cohesive energy', MINUS + '1.549 eV/atom', MINUS + '1.51 eV/atom (experiment) [@kittel]'],
    ['Basal GSFE, maximum along the straight ⟨112̅0⟩ path', '320.7 mJ/m²',
     'Not comparable with the DFT minimum-energy path'],
    ['Basal intrinsic stacking-fault energy (partial path)', '30.0 mJ/m²', '37.6 mJ/m² (DFT, I₂ fault) [@dong]'],
    ['Basal unstable stacking-fault energy (partial path)', '125.6 mJ/m²', '89.6 mJ/m² (DFT, I₂ fault) [@dong]'],
    ['Liquid mixing enthalpy, Mg₆₀Zn₃₅Ca₅, 1000 K', MINUS + '7.5 ± 0.05 kJ/mol', 'No reference value located'],
    ['C11', '62.81 GPa', '≈ 63 GPa (NIST) [@becker; @hale]'],
    ['C33', '69.62 GPa', '≈ 70 GPa (NIST) [@becker; @hale]'],
    ['C12', '25.97 GPa', '≈ 25 GPa (NIST) [@becker; @hale]'],
    ['C13', '21.18 GPa', '≈ 21 GPa (NIST) [@becker; @hale]'],
    ['C44', '17.14 GPa', '≈ 17 GPa (NIST) [@becker; @hale]'],
    ['C66', '18.42 GPa', '≈ 19 GPa (NIST) [@becker; @hale]'],
    ['Young’s modulus, Hill average (polycrystal)', '48.3 GPa', 'Derived from the constants above [@hill]'],
]
T2 = build_table(doc, ['Property', 'This work', 'Literature / reference'], t2_rows, [3500, 1900, 3960],
                 font_pt=9.5, align_first_left=True)

# ---------------------------------------------------------------- front matter
put('title', FRONT.TITLE, AL.CENTER)
put_runs('auth', FRONT.AUTHORS, AL.CENTER)
put_runs('body', FRONT.AFFIL, AL.CENTER, size=10)
put_runs('body', FRONT.CORR, AL.CENTER, size=10)
put('H1', 'Abstract')
put('body', FRONT.ABSTRACT, AL.JUSTIFY)
put_runs('body', [('Keywords: ', 'b'), (FRONT.KEYWORDS, False)], AL.LEFT)

for mod in (t15_a, t15_b, t15_d):
    for kind, payload in mod.BLOCKS:
        if kind in ('H1', 'H2'):
            keep_next(put(kind, payload))
        elif kind == 'P':
            put('body', payload, AL.JUSTIFY)
        elif kind == 'EQ':
            _append(copy.deepcopy(EQ[payload]))
        elif kind == 'FIG':
            put_fig(*payload)
        elif kind == 'TBL':
            src, idx, capt = payload
            keep_next(put_caption(capt))
            _append(T2 if src == 'new:T2' else T14[idx])
            _append(OxmlElement('w:p'))
for kind, txt in t15_d.DECL:
    if kind == 'H1':
        keep_next(put('H1', txt))
    else:
        put('body', txt, AL.JUSTIFY)

# ---------------------------------------------------------------- citations -> numbers, CMS style [3,6] and [7–9]
order = []
for t in body.iter(qn('w:t')):
    if t.text and '[@' in t.text:
        for m in re.finditer(r'\[@([^\]]+)\]', t.text):
            for k in re.findall(r'@([a-z0-9_]+)', m.group(0)):
                if k not in REFS:
                    raise SystemExit('UNKNOWN KEY ' + k)
                if k not in order:
                    order.append(k)
NUM = {k: i + 1 for i, k in enumerate(order)}

def keys_to_num(m):
    ns = sorted(NUM[k] for k in re.findall(r'@([a-z0-9_]+)', m.group(0)))
    out, i = [], 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
            j += 1
        out += ['%d–%d' % (ns[i], ns[j])] if j - i >= 2 else [str(x) for x in ns[i:j + 1]]
        i = j + 1
    return '[' + ','.join(out) + ']'

for t in body.iter(qn('w:t')):
    if t.text and '[@' in t.text:
        t.text = re.sub(r'\[@[^\]]+\]', keys_to_num, t.text)

ref_head = put('H1', 'References')
for k in order:
    put('body', '[%d] %s' % (NUM[k], REFS[k]), AL.LEFT)

# ---------------------------------------------------------------- tables: wording and CMS styling
CELL_REPL = [('Yield strength σy (GPa)', 'Peak strength σmax (GPa)'), ('Δσy vs glass (%)', 'Δσmax vs glass (%)'),
             ('σy (GPa)', 'σmax (GPa)'), ('Yield strain (%)', 'Peak strain (%)'), ('εy (%)', 'εmax (%)'),
             ('(yield)', '(peak)')]
CELL_EXACT = {'vf_1': 'B5 (vf_1)', 'vf_3': 'B15 (vf_3)', 'vf_1 (5% VF)': 'B5 (5% VF)',
              'vf_2 (10% VF)': 'B10 (10% VF)', 'vf_3 (15% VF)': 'B15 (15% VF)'}

def restyle_tables(document):
    n = 0
    for tbl in document.tables:
        seen = set()
        for row in tbl.rows:
            for c in row.cells:
                if c._tc in seen:
                    continue
                seen.add(c._tc)
                s = c.text; new = CELL_EXACT.get(s.strip(), s)
                for a, b in CELL_REPL:
                    new = new.replace(a, b)
                if new != s:
                    set_cell(c, new); n += 1
        t = tbl._tbl
        for el in list(t.iter(qn('w:shd'))) + list(t.iter(qn('w:tcBorders'))):
            el.getparent().remove(el)
        tblPr = t.tblPr
        old = tblPr.find(qn('w:tblBorders'))
        if old is not None:
            tblPr.remove(old)
        tb = OxmlElement('w:tblBorders')
        for side in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
            e = OxmlElement('w:' + side)
            if side in ('top', 'bottom'):
                e.set(qn('w:val'), 'single'); e.set(qn('w:sz'), '8'); e.set(qn('w:space'), '0'); e.set(qn('w:color'), '000000')
            else:
                e.set(qn('w:val'), 'nil')
            tb.append(e)
        # tblBorders must precede tblLayout/tblCellMar/tblLook in the schema order
        anchor = None
        for tag in ('w:tblLayout', 'w:tblCellMar', 'w:tblLook', 'w:tblCaption', 'w:tblDescription'):
            anchor = tblPr.find(qn(tag))
            if anchor is not None:
                break
        if anchor is not None:
            anchor.addprevious(tb)
        else:
            tblPr.append(tb)
        hdr = t.findall(qn('w:tr'))[0]
        for tc in hdr.findall(qn('w:tc')):
            tcPr = tc.find(qn('w:tcPr'))
            if tcPr is None:
                tcPr = OxmlElement('w:tcPr'); tc.insert(0, tcPr)
            b = OxmlElement('w:tcBorders'); bo = OxmlElement('w:bottom')
            bo.set(qn('w:val'), 'single'); bo.set(qn('w:sz'), '4'); bo.set(qn('w:space'), '0'); bo.set(qn('w:color'), '000000')
            b.append(bo)
            w = tcPr.find(qn('w:tcW'))
            (w.addnext(b) if w is not None else tcPr.insert(0, b))
    return n

print('table cells reworded:', restyle_tables(doc))

# ---------------------------------------------------------------- sub/superscript pass
MARK = re.compile(r'_\{([^}]*)\}|\^\{([^}]*)\}')
IMPLICIT = [
    re.compile(r'(?<![A-Za-z])(x)(Zn|Ca|Mg|B)(?![a-z])'),
    re.compile(r'(σ)(max|15|vM)'),
    re.compile(r'(ε)(max)'),
    re.compile(r'(?<![A-Za-z0-9])(q)(4|6|l)(?![0-9A-Za-z])'),
    re.compile(r'(?<![A-Za-z])(T)(g|f|0)(?![A-Za-z0-9])'),
    re.compile(r'(ΔH)(mix)'),
    re.compile(r'(V)_(Ca|Zn|vor)'),
    re.compile(r'(α)(CaCa|ZnZn|AB|AA)'),
    re.compile(r'(?<![A-Za-z])(P)(AB)(?![A-Za-z])'),
    re.compile(r'(D²)(min)'),
    re.compile(r'(?<![A-Za-z])(C)(11|33|12|13|44|66)(?![0-9])'),
]

def segments(text):
    segs, pos = [], 0
    for m in MARK.finditer(text):
        segs.append((text[pos:m.start()], None))
        segs.append((m.group(1), 'subscript') if m.group(1) is not None else (m.group(2), 'superscript'))
        pos = m.end()
    segs.append((text[pos:], None))
    for pat in IMPLICIT:
        out = []
        for s, va in segs:
            if va is not None or not s:
                out.append((s, va)); continue
            p0 = 0
            for m in pat.finditer(s):
                out.append((s[p0:m.start()], None)); out.append((m.group(1), None)); out.append((m.group(2), 'subscript'))
                p0 = m.end()
            out.append((s[p0:], None))
        segs = out
    return [(s, va) for s, va in segs if s]

def rich_paragraph(p_el):
    for r in list(p_el.iter(qn('w:r'))):
        ts = r.findall(qn('w:t'))
        if len(ts) != 1 or not ts[0].text:
            continue
        segs = segments(ts[0].text)
        if len(segs) == 1 and segs[0][1] is None:
            continue
        for s, va in segs:
            nr = copy.deepcopy(r)
            nt = nr.find(qn('w:t')); nt.text = s; nt.set(XMLSP, 'preserve')
            if va:
                rpr = nr.find(qn('w:rPr'))
                if rpr is None:
                    rpr = OxmlElement('w:rPr'); nr.insert(0, rpr)
                e = rpr.find(qn('w:vertAlign'))
                if e is None:
                    e = OxmlElement('w:vertAlign'); rpr.append(e)
                e.set(qn('w:val'), va)
            r.addprevious(nr)
        r.getparent().remove(r)

def rich_document(document, stop_el=None):
    for el in list(document.element.body):
        if stop_el is not None and el is stop_el:
            break
        for p_el in ([el] if el.tag == qn('w:p') else list(el.iter(qn('w:p')))):
            rich_paragraph(p_el)

rich_document(doc, ref_head._p)

# ---------------------------------------------------------------- continuous line numbers for review
sect = body.find(qn('w:sectPr'))
ln = OxmlElement('w:lnNumType'); ln.set(qn('w:countBy'), '1'); ln.set(qn('w:restart'), 'continuous')
pgsz = sect.find(qn('w:pgMar'))
(pgsz.addnext(ln) if pgsz is not None else sect.append(ln))

doc.save(OUT)
json.dump({'order': order, 'num': NUM}, open('/tmp/w/b15/ref_order_15.json', 'w'))
print('saved', OUT, '| references:', len(order))
print('unused keys:', [k for k in REFS if k not in NUM])
