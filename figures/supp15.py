# -*- coding: utf-8 -*-
"""Draft 15 supplement: draft 14 supplement with the new title/authors, peak-strength wording, B5/B15 names,
updated section numbers, CMS table style and sub/superscripts."""
import sys, io, contextlib
sys.path.insert(0, '/tmp/w/b15'); sys.path.insert(1, '/tmp/w/b14')
with contextlib.redirect_stdout(io.StringIO()):
    import build15 as B
import docx
from lib import settext
import t15_front as FRONT

d = docx.Document('/mnt/user-data/outputs/draft_14_supplement.docx')
P = d.paragraphs
settext(P[1], FRONT.TITLE)
settext(P[2], 'Ridhin V Raj, Aparna Zagabathuni (Department of Materials Science and Engineering, National Institute '
              'of Technology Calicut, India)')

TXT = [
    ('yield strain', 'peak strain'), ('(yield)', '(peak)'), ('Section 2.8)', 'Section 2.12)'),
    ('(59 in vf_1 and 177 in vf_3)', '(59 in B5, case vf_1, and 177 in B15, case vf_3)'),
    ('(1.20 in vf_1 and 1.14 in vf_3 at 7% strain)', '(1.20 in B5 and 1.14 in B15 at 7% strain)'),
    ('Supplementary Table S1. Mechanical properties of all 48 simulations.',
     'Supplementary Table S1. Mechanical properties of all 48 simulations. σmax is the peak strength and εmax the '
     'peak strain (Section 2.5 of the main text). B5 and B15 of the main text are cases vf_1 and vf_3.'),
]
n = 0
for p in P:
    s = p.text; new = s
    for a, b in TXT:
        new = new.replace(a, b)
    if new != s:
        settext(p, new); n += 1
print('supplement paragraphs edited:', n)
print('supplement table cells reworded:', B.restyle_tables(d))
B.rich_document(d)
left = [p.text[:80] for p in d.paragraphs if 'ield' in p.text]
print('remaining "yield":', left)
d.save('/tmp/w/b15/draft_15_supplement.docx')
print('saved')
