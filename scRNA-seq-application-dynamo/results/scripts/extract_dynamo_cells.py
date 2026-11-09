#!/usr/bin/env python3
"""extract_dynamo_cells.py -- turn the tutorial's nbconvert HTML into a runnable plain-Python
driver, one code cell at a time, in notebook order.

Why: the tutorial ships as scRNAseq-application_dynamo.html (nbconvert export of the teacher's
notebook). To run it "as written" on a login node we need the 34 code cells as a script; taking
them straight out of the HTML avoids re-typing and keeps the diff against the tutorial auditable.

The emitted file contains ONLY tutorial code plus a 6-line figure-flush harness after each cell
(the notebook shows figures inline; a script must save them). Every deviation from the notebook
text is applied afterwards by hand and is listed in the driver's @@NOTE header.

usage: python extract_dynamo_cells.py <tutorial.html> <out.py>
"""
import re
import sys
from html.parser import HTMLParser

VOID = {'br', 'img', 'hr', 'meta', 'link', 'input', 'source', 'col'}


class CodeCells(HTMLParser):
    """Collect the source of every <div class="jp-Cell jp-CodeCell"> input <pre> block."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.cells = []          # [(cell_index, source_text)]
        self.idx = -1
        self.buf = None

    def handle_starttag(self, tag, attrs):
        if tag in VOID:
            return
        cls = (dict(attrs).get('class') or '').split()
        info = {'tag': tag, 'cls': cls}
        if tag == 'div' and 'jp-Cell' in cls and 'jp-CodeCell' in cls:
            self.idx += 1
            self.cells.append((self.idx, []))
            info['incell'] = self.cells[-1]
        if tag == 'div' and 'jp-InputArea' in cls:
            info['input'] = True
        if tag == 'div' and 'highlight' in cls:
            info['hl'] = True
        if tag == 'pre':
            incell = self.top('incell')
            if incell is not None and self.top('hl') and self.top('input'):
                self.buf = incell[1]
                info['buf'] = True
        self.stack.append((tag, info))

    def top(self, key):
        for _t, info in reversed(self.stack):
            if key in info:
                return info[key]
        return None

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        while self.stack:
            t, info = self.stack.pop()
            if tag == 'pre' and info.get('buf'):
                self.buf = None
            if t == tag:
                break

    def handle_data(self, data):
        if self.buf is not None:
            self.buf.append(data)


if __name__ == '__main__':
    html_path, out_path = sys.argv[1], sys.argv[2]
    p = CodeCells()
    p.feed(open(html_path, encoding='utf-8', errors='replace').read())
    cells = [(i, re.sub(r'\n+$', '', ''.join(parts))) for i, parts in p.cells]
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('#!/usr/bin/env python\n')
        f.write('# AUTO-EXTRACTED from the tutorial HTML by extract_dynamo_cells.py -- do not hand-edit\n')
        f.write('# the cell bodies; apply deviations as a documented patch at the top instead.\n')
        f.write('# @@CELLS %d code cells\n' % len(cells))
        f.write('\n')
        f.write('# ===== harness: flush any figure left open by the cell that just ran =====\n')
        f.write('import os\n')
        f.write('import matplotlib\n')
        f.write('matplotlib.use("Agg")\n')
        f.write('import matplotlib.pyplot as _plt\n')
        f.write('_FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figure")\n')
        f.write('os.makedirs(_FIGDIR, exist_ok=True)\n')
        f.write('_INLINE = []\n')
        f.write('def _flush(cell, tag):\n')
        f.write('    for num in _plt.get_fignums():\n')
        f.write('        fig = _plt.figure(num)\n')
        f.write('        if not fig.get_axes():\n')
        f.write('            continue\n')
        f.write('        name = "inline_cell%02d_%s_%02d.png" % (cell, tag, num)\n')
        f.write('        fig.savefig(os.path.join(_FIGDIR, name), dpi=200, bbox_inches="tight")\n')
        f.write('        _INLINE.append(name)\n')
        f.write('        print("[harness] wrote figure/%s" % name)\n')
        f.write('    _plt.close("all")\n')
        f.write('\n')
        for idx, src in cells:
            f.write('\n# ===================== cell %d (verbatim) =====================\n' % idx)
            f.write(src.rstrip('\n') + '\n')
            f.write('_flush(%d, "auto")\n' % idx)
        f.write('\n# ===================== end of run =====================\n')
        f.write('print("[harness] inline figures captured: %d" % len(_INLINE))\n')
        f.write('for _n in _INLINE:\n')
        f.write('    print("[harness]   " + _n)\n')
    print('%s: %d code cells -> %s' % (html_path, len(cells), out_path))
