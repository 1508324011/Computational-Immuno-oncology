#!/usr/bin/env python3
"""extract_html_images.py -- pull the base64 `data:image/png` payloads out of an nbconvert HTML export.

nbconvert embeds every cell output as a data URI, so the only way to see what the tutorial's figures
looked like (and to compare against them) is to decode them.  Images are written in document order.

usage: python extract_html_images.py <tutorial.html> <out_dir>
"""
import base64
import hashlib
import os
import re
import sys

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
html = open(src, encoding='utf-8', errors='replace').read()
pat = re.compile(r'data:image/(png|jpeg|jpg|svg\+xml);base64,([A-Za-z0-9+/=\s]+)')
seen, rows = set(), []
for i, m in enumerate(pat.finditer(html), 1):
    ext = 'png' if m.group(1).startswith('png') else ('jpg' if 'jp' in m.group(1) else 'svg')
    blob = base64.b64decode(re.sub(r'\s+', '', m.group(2)))
    h = hashlib.md5(blob).hexdigest()
    dup = h in seen
    seen.add(h)
    name = '%02d_%s.%s' % (i, h[:8], ext) if not dup else '%02d_%s_dup.%s' % (i, h[:8], ext)
    open(os.path.join(out, name), 'wb').write(blob)
    rows.append((name, len(blob), dup))
    print('%-26s %8d bytes%s' % (name, len(blob), '  (duplicate payload)' if dup else ''))
print('\n%d images, %d unique payloads -> %s' % (len(rows), len(seen), out))
