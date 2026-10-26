#!/usr/bin/env python
"""run_cells.py -- execute the tutorial's Python cells in ONE persistent namespace,
notebook-style: a failing cell is recorded and the next cell still runs.
Nothing is auto-corrected.

usage: python run_cells.py <cells_file> <out_dir> <label>
"""
import io
import re
import sys
import time
import traceback
import warnings
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

cells_file, out_dir, label = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
out_dir.mkdir(parents=True, exist_ok=True)
logf = open(out_dir / f"{label}_cell_log.txt", "w", buffering=1)


def say(*parts):
    txt = " ".join(str(p) for p in parts)
    print(txt, flush=True)
    logf.write(txt + "\n")


src = Path(cells_file).read_text().splitlines()
idx = [(i, int(m.group(1))) for i, m in
       ((i, re.match(r"^#@@CELL (\d+)$", ln)) for i, ln in enumerate(src)) if m]
if not idx:
    raise SystemExit("no '#@@CELL <id>' sentinels found")
cells = []
for k, (i, cid) in enumerate(idx):
    end = idx[k + 1][0] if k + 1 < len(idx) else len(src)
    cells.append((cid, "\n".join(src[i + 1:end])))

import scanpy as sc
import anndata as ad

say("### driver: run_cells.py | label=", label)
say("### python:", sys.version.split()[0], "| scanpy:", sc.__version__,
    "| anndata:", ad.__version__, "| cells:", len(cells))
say("### cwd:", Path.cwd())
say("### start:", time.strftime("%Y-%m-%d %H:%M:%S"))

ns = {"__name__": "__main__"}
rows = []


def probe():
    out = []
    adata = ns.get("adata")
    if adata is not None:
        out.append(f"adata: {adata.n_obs} cells x {adata.n_vars} genes")
        if "raw" in adata.__dict__ and adata.raw is not None:
            out.append(f"adata.raw: {adata.raw.n_vars} genes")
        for col in ("leiden", "n_genes_by_counts", "pct_counts_mt"):
            if col in adata.obs:
                s = adata.obs[col]
                if s.dtype.name == "category":
                    out.append(f"obs['{col}']: {len(s.cat.categories)} categories "
                               f"[{'; '.join(list(s.cat.categories)[:10])}]")
                else:
                    out.append(f"obs['{col}']: median={s.median():.3f} max={s.max():.1f}")
        for key in ("pca", "neighbors", "umap"):
            if key in adata.obsm or key in adata.uns:
                pass
        if "X_pca" in adata.obsm:
            out.append(f"obsm['X_pca']: {adata.obsm['X_pca'].shape}")
        if "X_umap" in adata.obsm:
            out.append(f"obsm['X_umap']: {adata.obsm['X_umap'].shape}")
        if "leiden" in adata.uns.get("rank_genes_groups", {}).get("params", {}):
            pass
        if "rank_genes_groups" in adata.uns:
            out.append("uns['rank_genes_groups']: present "
                       f"(method={adata.uns['rank_genes_groups'].get('params', {}).get('method')})")
    if "adata.raw" in ns and ns["adata.raw"] is not None:
        out.append(f"adata.raw: {ns['adata.raw'].n_vars} genes")
    return out


for cid, code in cells:
    say("")
    say(f"===== CELL {cid} =====")
    for ln in code.splitlines():
        say("  |", ln)
    buf_out, buf_err = io.StringIO(), io.StringIO()
    warns, err, t0 = [], None, time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            with redirect_stdout(buf_out), redirect_stderr(buf_err):
                exec(compile(code, f"<cell {cid}>", "exec"), ns)
        except BaseException as exc:  # noqa: BLE001 - we want everything, incl. SyntaxError
            err = exc
        warns = [f"{w.category.__name__}: {w.message}" for w in caught]
    elapsed = time.time() - t0
    o, e = buf_out.getvalue().strip(), buf_err.getvalue().strip()
    if o:
        say("  --stdout--")
        for ln in o.splitlines():
            say("  >", ln)
    if e:
        say("  --stderr--")
        for ln in e.splitlines():
            say("  >", ln)
    uniq = list(dict.fromkeys(warns))
    if uniq:
        say(f"  --warnings({len(uniq)})--")
        for w in uniq:
            say("  !", w)
    if err is not None:
        say("  --ERROR--", f"{type(err).__name__}: {err}")
        tb = traceback.format_exception(type(err), err, err.__traceback__)
        for ln in "".join(tb).splitlines()[-6:]:
            say("  ~", ln)
    st = probe()
    if st:
        say("  --state--")
        for s in st:
            say("  #", s)
    say(f"  --status-- {'OK' if err is None else 'ERROR'} in {elapsed:.2f}s")
    rows.append(dict(cell=cid, status="OK" if err is None else "ERROR",
                     seconds=round(elapsed, 2), n_warn=len(uniq),
                     error="" if err is None else f"{type(err).__name__}: {err}",
                     warn_first=uniq[0] if uniq else "",
                     state=" ;; ".join(st)))

import csv
with open(out_dir / f"{label}_cell_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

say("")
say("===== SUMMARY =====")
for r in rows:
    say(f"cell {r['cell']:<3d} {r['status']:<5s} {r['seconds']:>7.2f}s "
        f"warn={r['n_warn']} {r['error'][:120]}")
say("### end:", time.strftime("%Y-%m-%d %H:%M:%S"))
logf.close()
