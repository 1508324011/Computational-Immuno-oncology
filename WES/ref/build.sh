#!/usr/bin/env bash
# Build WES_roadmap_corrected.{html,pdf} from the markdown source.
#
# Why the post-processing exists (2026-09-21):
#   1. pandoc wraps every table in `#figure(align(center)[#table(...)], kind: table)`.
#      typst figures are NOT breakable across pages: when a tall table lands near a
#      page bottom, typst 0.15.1 compresses the remaining rows to ~one line each and
#      the rows overlap (reproduced minimally: figure-wrap -> 135 overlaps, no-wrap -> 0).
#      Fix: strip the figure wrapper so tables break across pages normally.
#   2. Inline `code` tokens >=12 chars (paths, flags, filenames) need break
#      opportunities: typst breaks raw spans after '/' but NOT after '_' '.' '-'.
#      Verified empirically (typst 0.15.1): a real U+200B zero-width space inside a
#      raw span IS honored as a break point. (Do NOT use a re replacement template
#      like r'\1\u200b' -- Python 3.6 silently inserts the literal text "\u200b";
#      build with chr(0x200b) instead.) Fenced code blocks are left untouched.
#   3. Narrow equal-percentage columns overflow even after breaking, so tables get
#      8pt cell text and tighter insets via a pandoc -H header include.
#
# Verification (automatic): pdftotext -bbox + word-bbox overlap / margin-overflow
# scan over every page. Known-benign artifacts are filtered (hanging punctuation on
# justified lines, fullwidth-quote kerning). Build fails loudly if anything remains.
set -euo pipefail
cd "$(dirname "$0")"

MD=WES_roadmap_corrected.md
source ~/miniconda3/etc/profile.d/conda.sh
conda activate docs   # pandoc 3.11, typst 0.15.1, python 3.6 (system-compatible only!)

# ---------------------------------------------------------------- HTML -------
pandoc "$MD" -o WES_roadmap_corrected.html --standalone --toc --toc-depth=2 \
  --metadata title="WES 全流程分析路线图 · 矫正增强版 (Corrected & Annotated Edition)" \
  --metadata lang=zh-CN

# ---------------------------------------------------------------- PDF --------
# NOTE: no --metadata lang for the PDF on purpose: lang=en keeps typst hyphenation
# (and thus clean wrapping of long English words) enabled; CJK wraps regardless.
HDR=$(mktemp /tmp/corrected_hdr_XXXX.typ)
TYP=$(mktemp /tmp/corrected_XXXX.typ)
BBOX=$(mktemp /tmp/bbox_XXXX.html)
trap 'rm -f "$TYP" "$BBOX" "$HDR"' EXIT

cat > "$HDR" <<'TYPST'
// denser tables: 8pt cells + tighter insets (equal-% columns are narrow)
#set table(inset: (x: 4pt, y: 4pt))
#show table.cell: set text(size: 8pt)
TYPST

pandoc "$MD" -t typst -s -o "$TYP" \
  --toc --toc-depth=2 \
  --metadata title="WES 全流程分析路线图 · 矫正增强版" \
  -V mainfont="Noto Sans CJK SC" -V fontsize=9pt \
  -H "$HDR"

python3 - "$TYP" <<'PYEOF'
import re, sys
path = sys.argv[1]
src = open(path, encoding='utf-8').read()
lines = src.split('\n')
out, i, stripped = [], 0, 0
while i < len(lines):
    if lines[i] == '#figure(' and i + 1 < len(lines) and lines[i + 1].startswith('  align(center)[#table('):
        out.append(lines[i + 1].replace('  align(center)[', '', 1))  # -> '  #table('
        i += 2
        while i < len(lines):
            if lines[i] == '  )]':
                out.append('  )')
                i += 1
                if i < len(lines) and lines[i] == '  , kind: table':
                    i += 1
                if i < len(lines) and lines[i] == ')':
                    i += 1
                stripped += 1
                break
            out.append(lines[i]); i += 1
        continue
    out.append(lines[i]); i += 1
src = '\n'.join(out)
sys.stderr.write('figures unwrapped: %d\n' % stripped)

ZWSP = chr(0x200b)  # real zero-width space (NOT the literal text '\u200b'!)
parts = re.split(r'(?m)(^```.*?^```$)', src, flags=re.S)
for k in range(0, len(parts), 2):  # even = prose/tables, odd = fenced code blocks
    p = parts[k]
    # (a) plain tokens: pandoc escapes '_' as '\_' in typst markup -- e.g. table
    # headers like PERCENT\_DUPLICATION / MEAN\_TARGET\_COV. typst never breaks
    # after an escaped underscore, so give it a ZWSP break point.
    p = re.sub(r'((?:[A-Za-z0-9]|\\_){11,})',
               lambda m: re.sub(r'(\\_)', lambda mm: mm.group(1) + ZWSP, m.group(1)),
               p)
    # (b) inline raw (backtick) tokens >=12 chars: typst breaks raws after '/'
    # but not after '_' '.' '-' -- insert ZWSPs after those separators.
    p = re.sub(r'`([^`\n]{12,})`',
               lambda m: '`' + re.sub(r'([/_.>\-])', lambda mm: mm.group(1) + ZWSP, m.group(1)) + '`',
               p)
    parts[k] = p
open(path, 'w', encoding='utf-8').write(''.join(parts))
sys.stderr.write('ZWSP break points inserted\n')
PYEOF

typst compile "$TYP" WES_roadmap_corrected.pdf
echo "== typst compile OK =="

# --------------------------------------------------------- verify: overlaps --
pdftotext -bbox WES_roadmap_corrected.pdf "$BBOX"
python3 - "$BBOX" <<'PYEOF'
import re, sys
html = open(sys.argv[1], encoding='utf-8').read()
pages = re.findall(r'<page width="([\d.]+)" height="([\d.]+)">(.*?)</page>', html, re.S)
total, report = 0, []
for pi, (pw, ph, body) in enumerate(pages, 1):
    words = [(float(a), float(b), float(c), float(d), t.replace(chr(0x200b), '~')) for a, b, c, d, t in
             re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>', body)]
    words.sort(key=lambda w: (round(w[1], 1), w[0]))
    n, bad, over = len(words), [], []
    for x in range(n):
        x0, y0, x1, y1, t = words[x]
        if x1 > 528 and re.search(r'[0-9A-Za-z]', t[-1:]):  # >6pt past the 522pt margin
            over.append('%s[x1=%.0f]' % (t[:24], x1))
        for y in range(x + 1, n):
            a0, b0, a1, b1, u = words[y]
            if b0 >= y1: break
            ox = min(x1, a1) - max(x0, a0); oy = min(y1, b1) - max(y0, b0)
            # real overlap: >2.5pt wide, >45% of the smaller glyph height, >15% of the
            # smaller box area. The 2.5pt floor filters fullwidth-quote/kerning artifacts.
            if ox > 2.5 and oy > 0.45 * min(y1 - y0, b1 - b0):
                frac = (ox * oy) / max(1e-9, min((x1 - x0) * (y1 - y0), (a1 - a0) * (b1 - b0)))
                # CJK fullwidth punctuation: the glyph box abuts the preceding text but
                # its ink sits in the right half -- a designed overlap, not a defect.
                cjk_punct = set('（）。，：；！？、。「」『』“”‘’—…·')
                if t in cjk_punct or u in cjk_punct:
                    continue
                if frac > 0.15: bad.append('%.2f %r<->%r' % (frac, t[:24], u[:24]))
            if b0 - y0 > 60: break
    if bad or over:
        report.append('page %d: overlaps=%d %s | right-overflow=%d %s'
                      % (pi, len(bad), bad[:3], len(over), over[:3]))
        total += len(bad) + len(over)
print('pages: %d, problem pages: %d, total flags: %d' % (len(pages), len(report), total))
for r in report: print(' ', r)
sys.exit(1 if total else 0)
PYEOF
echo "== overlap scan clean =="

# ---------------------------------------------------------------- fonts ------
echo "== fonts used =="
pdffonts WES_roadmap_corrected.pdf 2>/dev/null | tail -n +3 | awk '{print "  ", $1}' | sed 's/^[^ ]*//' | sort -u | head -8
echo "== build done: $(pdfinfo WES_roadmap_corrected.pdf 2>/dev/null | grep -E '^Pages') =="
