#!/usr/bin/env python
"""compare_figs.py -- quantitative, image-based comparison of our re-run figures
against the 28 reference PNGs that ship with the course (the teacher's 2025 output,
already renamed into class/scRNA-seq/*.png).

Per figure: presence, pixel size, aspect, 64x64 structural correlation, and -- for the
panel types where it is meaningful -- the number of distinct painted colours and their
pixel shares (a proxy for cell-type composition readable straight off a PNG).

usage: python compare_figs.py <run_root> <ref_dir> <out_tsv> [seurat_subdir] [py_subdir]
        seurat_subdir default: rr_fixed (the corrected run; pass rr for the as-written run)
        py_subdir     default: py_fixed
"""
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image

SEURAT = {
    "violin_mt.png": ("violin_mt.png", "violin"),
    "VariableFeaturePlot.png": ("VariableFeaturePlot.png", "scatter"),
    "PCA_VizDimLoadings.png": ("PCA_VizDimLoadings.png", "bar"),
    "PCA.png": ("PCA.png", "line"),
    "PCA_dimplot.png": ("PCA_dimplot.png", "umap"),
    "PCA_heatmap_ALL.png": ("PCA_heatmap_ALL.png", "heatmap"),
    "UMAP.png": ("UMAP.png", "umap"),
    "marker_Features.png": ("marker_Features.png", "featuregrid"),
    "Marker_Heatmap.png": ("Marker_Heatmap.png", "heatmap"),
    "marker_Features_2.png": ("marker_Features_3.png", "featuregrid"),
    "cell annotations.png": ("cell annotations.png", "umap"),
    "Naive_Memory_CD4_DEG.png": ("Naive_Memory_CD4_DEG.png", "volcano"),
    "UMAP_without_harmony.png": ("UMAP_without_harmony.png", "umap"),
    "UMAP_with_harmony.png": ("UMAP_with_harmony.png", "umap"),
}
SCANPY = {
    "figures/highest_expr_genesscanpy_highest_expr_genes.png": ("scanpy_highest_expr_genes.png", "bar"),
    "figures/violinscanpy_violin_plot_QC.png": ("scanpy_violin_plot_QC.png", "violin"),
    "figures/scatterscanpy_scatter_pct_counts_mt.png": ("scanpy_scatter_pct_counts_mt.png", "scatter"),
    "figures/scatterscanpy_scatter_n_genes_by_counts.png":
        ("scanpy_scatter_n_genes_by_counts.png", "scatter"),
    "figures/filter_genes_dispersionscanpy_HVG.png": ("scanpy_HVG.png", "scatter"),
    "figures/pcascanpy_PCA.png": ("scanpy_PCA.png", "umap"),
    "figures/pca_variance_ratioscanpy_PCA_ratio.png": (None, "line"),
    "figures/umapscanpy_markers.png": ("scanpy_markers.png", "featuregrid"),
    "figures/rank_genes_groups_leidenscanpy_rank_genes_groups.png":
        ("scanpy_rank_genes_groups.png", "rankplot"),
    "figures/rank_genes_groups_leiden_0scanpy_rank_genes_groups_violin.png":
        ("scanpy_rank_genes_groups_violin.png", "violin"),
    "figures/dotplot_scanpy_rank_genes_groups_dotplot.png":
        ("scanpy_rank_genes_groups_dotplot.png", "dotplot"),
    "figures/umapscanpy_final_umap.png": ("scanpy_final_umap.png", "umap"),
}


def structural(path, size=64):
    im = Image.open(path).convert("L").resize((size, size), Image.BILINEAR)
    a = np.asarray(im, dtype=np.float64)
    a -= a.mean()
    n = np.linalg.norm(a)
    return a / n if n else a


def colours(path, min_share=0.004, sat=40, val=40):
    im = Image.open(path).convert("RGB")
    im.thumbnail((900, 900))
    a = np.asarray(im, dtype=np.int16)
    mx, mn = a.max(2), a.min(2)
    painted = (mx - mn > sat) & (mx > val)
    if painted.sum() == 0:
        return 0, {}
    keys, counts = np.unique((a[painted] // 24) * 24, axis=0, return_counts=True)
    tot = counts.sum()
    out = {}
    for i in np.argsort(-counts):
        share = counts[i] / tot
        if share < min_share:
            break
        out["#%02x%02x%02x" % tuple(keys[i])] = round(float(share), 3)
    return len(out), out


def main():
    run_root, ref_dir, out_tsv = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    seurat_sub = sys.argv[4] if len(sys.argv) > 4 else "rr_fixed"
    py_sub = sys.argv[5] if len(sys.argv) > 5 else "py_fixed"
    rows = []
    for base, table in ((seurat_sub, SEURAT), (py_sub, SCANPY)):
        for rel, (ref, kind) in table.items():
            op = run_root / base / rel
            rp = ref_dir / ref if ref else None
            row = {"run": base, "ours": rel, "reference": ref or "(not produced by the tutorial)",
                   "kind": kind, "ours_exists": op.exists(),
                   "ref_exists": bool(rp and rp.exists())}
            if op.exists():
                io_ = Image.open(op)
                row["ours_px"] = f"{io_.size[0]}x{io_.size[1]}"
                row["ours_kb"] = round(op.stat().st_size / 1024, 1)
            if rp and rp.exists():
                ir = Image.open(rp)
                row["ref_px"] = f"{ir.size[0]}x{ir.size[1]}"
                row["ref_kb"] = round(rp.stat().st_size / 1024, 1)
            if op.exists() and rp and rp.exists():
                row["struct_corr"] = round(float((structural(op) * structural(rp)).sum()), 4)
                if kind in ("umap", "featuregrid", "scatter", "volcano", "dotplot", "rankplot"):
                    no, co = colours(op)
                    nr, cr = colours(rp)
                    row["ours_colours"] = no
                    row["ref_colours"] = nr
                    row["ours_shares"] = ";".join(f"{k}:{v}" for k, v in co.items())
                    row["ref_shares"] = ";".join(f"{k}:{v}" for k, v in cr.items())
            rows.append(row)

    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['run']:10s} {r['ours'][-56:]:56s} ours={str(r['ours_exists'])[0]} "
              f"ref={str(r['ref_exists'])[0]} {r.get('ours_px','-'):>11s} vs "
              f"{r.get('ref_px','-'):>11s} corr={r.get('struct_corr','-')}")


if __name__ == "__main__":
    main()
