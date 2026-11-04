#!/usr/bin/env python
"""compare_figs_velo.py -- quantitative comparison of our scVelo re-run figures against the
11 reference PNGs that ship with this module (class/scRNA-seq-scVelo/*.png, i.e. the teacher's
own notebook output extracted from scRNAseq-scVelo.html).

Per figure: presence, pixel size, aspect ratio, 64x64 structural correlation (Pearson on the
mean-centred grayscale thumbnails -- robust to size/dpi differences, sensitive to layout), and
the number of distinct painted colours plus their pixel shares (a proxy for how many cell
clusters are painted and in what proportion, readable straight off a PNG).

usage: python compare_figs_velo.py <run_root> <ref_dir> <out_tsv> [figs_subdir]
        figs_subdir default: figures
"""
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image

FIGURES = [
    ("01_proportions.png", "cell 50  scv.pl.proportions -- spliced/unspliced share", "pie+bar"),
    ("02_velocity_embedding_stream.png", "cell 63  scv.pl.velocity_embedding_stream", "stream"),
    ("03_velocity_embedding.png", "cell 66  scv.pl.velocity_embedding (arrows)", "quiver"),
    ("04_phaseportrait_markers.png", "cell 69  scv.pl.velocity(Cpe,Gnao1,Ins2,Adk)", "phase"),
    ("05_cpe_groups_velocity.png", "cell 71  scv.pl.scatter(Cpe, clusters+velocity)", "featuregrid"),
    ("06_top5_ngn3_high_ep.png", "cell 75  scv.pl.scatter(top5 Ngn3 high EP)", "featuregrid"),
    ("07_top5_pre_endocrine.png", "cell 75  scv.pl.scatter(top5 Pre-endocrine)", "featuregrid"),
    ("08_cellcycle_gradients.png", "cell 79  scv.pl.scatter(color_gradients S/G2M)", "gradient"),
    ("09_phase_marker_genes.png", "cell 81  scv.pl.scatter(phase marker genes)", "featuregrid"),
    ("10_phaseportrait_hells_top2a.png", "cell 83  scv.pl.velocity(Hells,Top2a)", "phase"),
    ("11_velocity_embedding_stream_repeat.png", "cell 84  scv.pl.velocity_embedding_stream (repeat)", "stream"),
]

COLOUR_KINDS = {"stream", "quiver", "featuregrid", "gradient"}

# The corrected pass does not have one panel per tutorial plot: the duplicated stream figure is
# replaced by latent time / confidence, and the phase portraits exist twice (stochastic, as the
# tutorial draws them, and dynamical, with the fitted kinetics). This table says which reference
# PNG each corrected panel is comparable to; None = no counterpart in the teacher's export.
CORRECTED_PAIRS = [
    ("corrected_01_proportions.png", "01_proportions.png", "proportions", "pie+bar"),
    ("corrected_02_velocity_embedding_stream.png", "02_velocity_embedding_stream.png", "stream", "stream"),
    ("corrected_03_velocity_embedding.png", "03_velocity_embedding.png", "arrows", "quiver"),
    ("corrected_04_phaseportrait_markers_stochastic.png", "04_phaseportrait_markers.png", "phase portraits, stochastic", "phase"),
    ("corrected_05_phaseportrait_markers_dynamical.png", "04_phaseportrait_markers.png", "phase portraits, dynamical (with fits)", "phase"),
    ("corrected_06_cpe_groups_velocity.png", "05_cpe_groups_velocity.png", "Cpe groups+velocity", "featuregrid"),
    ("corrected_07_top5_ngn3_high_ep.png", "06_top5_ngn3_high_ep.png", "top5 Ngn3 high EP", "featuregrid"),
    ("corrected_08_top5_pre_endocrine.png", "07_top5_pre_endocrine.png", "top5 Pre-endocrine", "featuregrid"),
    ("corrected_09_cellcycle_gradients_original_scores.png", "08_cellcycle_gradients.png", "cell cycle, original scores", "gradient"),
    ("corrected_10_phase_marker_genes.png", "09_phase_marker_genes.png", "phase marker genes", "featuregrid"),
    ("corrected_11_phaseportrait_hells_top2a_dynamical.png", "10_phaseportrait_hells_top2a.png", "Hells/Top2a, dynamical", "phase"),
    ("corrected_12_latent_time_and_confidence.png", None, "latent time + confidence (replaces the duplicate)", "gradient"),
    ("corrected_13_velocity_stream_latent_time.png", "02_velocity_embedding_stream.png", "stream coloured by latent time", "stream"),
]


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


def ink_mask(path, size=256, white=243):
    """Binary 'something was drawn here' mask: any channel below the near-white level.

    For low-contrast panels (the cell-cycle gradient plot is mostly white background) the
    grayscale structural correlation is dominated by faint background differences, so the
    painted-area masks are compared as an additional, contrast-insensitive metric.
    """
    im = Image.open(path).convert("RGB").resize((size, size), Image.BILINEAR)
    a = np.asarray(im, dtype=np.uint8)
    return a.min(axis=2) < white


def mask_iou(pa, pb):
    ma, mb = ink_mask(pa), ink_mask(pb)
    inter = np.logical_and(ma, mb).sum()
    union = np.logical_or(ma, mb).sum()
    return (inter / union if union else 0.0), ma.mean(), mb.mean()


def main():
    run_root, ref_dir, out_tsv = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    sub = sys.argv[4] if len(sys.argv) > 4 else "figures"
    which = sys.argv[5] if len(sys.argv) > 5 else "aswritten"
    pairs = [(n, r, w, k) for (n, w, k) in FIGURES for r in [n]] if which == "aswritten" else CORRECTED_PAIRS
    rows = []
    for name, ref, what, kind in pairs:
        op = run_root / sub / name
        rp = (ref_dir / ref) if ref else None
        row = {
            "figure": name,
            "tutorial_call": what,
            "kind": kind,
            "ours_exists": op.exists(),
            "ref_exists": bool(rp and rp.exists()),
        }
        if op.exists():
            io_ = Image.open(op)
            row["ours_px"] = "%dx%d" % io_.size
            row["ours_kb"] = round(op.stat().st_size / 1024, 1)
        if rp and rp.exists():
            ir = Image.open(rp)
            row["ref_px"] = "%dx%d" % ir.size
            row["ref_kb"] = round(rp.stat().st_size / 1024, 1)
        if op.exists() and rp and rp.exists():
            io_, ir = Image.open(op), Image.open(rp)
            row["ours_aspect"] = round(io_.size[0] / io_.size[1], 3)
            row["ref_aspect"] = round(ir.size[0] / ir.size[1], 3)
            row["struct_corr"] = round(float((structural(op) * structural(rp)).sum()), 4)
            iou, sha, shr = mask_iou(op, rp)
            row["ink_iou"] = round(float(iou), 4)
            row["ours_ink_share"] = round(float(sha), 4)
            row["ref_ink_share"] = round(float(shr), 4)
            if kind in COLOUR_KINDS:
                no, co = colours(op)
                nr, cr = colours(rp)
                row["ours_colours"] = no
                row["ref_colours"] = nr
                row["ours_shares"] = ";".join("%s:%s" % (k, v) for k, v in co.items())
                row["ref_shares"] = ";".join("%s:%s" % (k, v) for k, v in cr.items())
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
        print(
            "%-40s ours=%s ref=%s %11s vs %11s  aspect %.2f/%.2f  corr=%s  ink_iou=%s"
            % (
                r["figure"],
                "y" if r["ours_exists"] else "n",
                "y" if r["ref_exists"] else "n",
                r.get("ours_px", "-"),
                r.get("ref_px", "-"),
                r.get("ours_aspect", 0),
                r.get("ref_aspect", 0),
                r.get("struct_corr", "-"),
                r.get("ink_iou", "-"),
            )
        )


if __name__ == "__main__":
    main()
