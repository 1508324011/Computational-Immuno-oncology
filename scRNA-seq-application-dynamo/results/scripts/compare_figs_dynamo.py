#!/usr/bin/env python
"""compare_figs_dynamo.py -- quantitative comparison of our Dynamo re-run figures against the
teacher's reference renders.

Two reference families are compared, as agreed for this module:
  * teacher2026 : /lustre1/share/scRNA-adv-2026/figure/dynamo/ (11 PNGs, produced 2026-11-08 by
                  the staff's scripted re-run of the same notebook, dpi 300, includes the three
                  Jacobian panels and the two GO panels that the tutorial itself never saves)
  * teacherhtml : the six figures embedded in the tutorial's own nbconvert export
                  (scRNAseq-application_dynamo.html -> ref/figures/teacherhtml/), i.e. exactly
                  what a student sees inline

Per pair we report presence, pixel size, aspect ratio, a 64x64 structural correlation (Pearson on
mean-centred grayscale thumbnails: robust to dpi differences, sensitive to layout/geometry), an
ink-mask IoU (contrast-insensitive, catches "same structures, different everything else"), and the
number of distinct painted colours with their pixel shares (a proxy for how many clusters are
painted and in what proportions, readable straight off a PNG).

usage: python compare_figs_dynamo.py <run_root> <ref_root> <out_tsv> [--figs SUBDIR]
      ref_root must contain teacher2026/ and teacherhtml/
"""
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image

# (ours, teacher2026 ref, teacherhtml ref, tutorial cell, what the call is, kind)
PAIRS = [
    ("01_velocity_streamline_dyn_savefig.png", "01_velocity_streamline_dyn_savefig.png", None,
     "cell 14", "dyn.pl.streamline_plot(clusters, umap)", "stream"),
    ("02_cell_wise_velocity_dyn_savefig.png", "02_cell_wise_velocity_dyn_savefig.png", None,
     "cell 15", "dyn.pl.cell_wise_vectors(clusters, umap)", "quiver"),
    ("03_combined_plots.png", "03_combined_plots.png", "cell030_01.png",
     "cell 18", "4-panel: ddhodge potential / speed / divergence / acceleration", "composite"),
    ("04_cell_state_dyn_savefig.png", "04_cell_state_dyn_savefig.png", None,
     "cell 19", "dyn.pl.topography(clusters)", "composite"),
    ("05_cell_state_transition_map_dyn_savefig.png", "05_cell_state_transition_map_dyn_savefig.png", None,
     "cell 20", "dyn.pd.state_graph + dyn.pl.state_graph", "graph"),
    ("06_reg_Ngn3_activate_Pax4_dyn_savefig.png", "06_reg_Ngn3_activate_Pax4_dyn_savefig.png", None,
     "cell 22", "dyn.pl.jacobian(Neurog3 -> Pax4)", "jacobian"),
    ("inline_cell23_auto_01.png", "07_reg_Arx_Pax4_dyn_savefig.png", "cell038_02.png",
     "cell 23", "dyn.pl.jacobian(Arx -> Pax4)  [tutorial does not save it]", "jacobian"),
    ("inline_cell23_auto_02.png", "08_reg_Pax4_Arx_dyn_savefig.png", "cell038_03.png",
     "cell 23", "dyn.pl.jacobian(Pax4 -> Arx)  [tutorial does not save it]", "jacobian"),
    ("inline_cell24_auto_01.png", "09_reg_Pax4_Ins2_dyn_savefig.png", "cell039_04.png",
     "cell 24", "dyn.pl.jacobian(Pax4 -> Ins2)  [tutorial does not save it]", "jacobian"),
    ("inline_cell32_auto_01.png", "10_Ductal_GO_enrichment_barplot.png", "cell048_05.png",
     "cell 32", "gseapy barplot(GO BP 2021, Ductal top100)", "bar"),
    ("inline_cell32_auto_02.png", "11_Ductal_GO_enrichment_dotplot.png", "cell048_06.png",
     "cell 32", "gseapy dotplot(same table)", "dot"),
]


def structural(path, size=64):
    """Mean-centred, unit-norm 64x64 grayscale thumbnail. None when the image carries no
    structure at all (a blank white panel), which some staff reference renders are."""
    im = Image.open(path).convert("L").resize((size, size), Image.BILINEAR)
    a = np.asarray(im, dtype=np.float64)
    a -= a.mean()
    n = np.linalg.norm(a)
    return (a / n) if n > 1e-9 else None


def ink_mask(path, size=256, white=243):
    im = Image.open(path).convert("RGB").resize((size, size), Image.BILINEAR)
    a = np.asarray(im, dtype=np.uint8)
    return a.min(axis=2) < white


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


def metrics(op, rp):
    row = {}
    oi, ri = Image.open(op), Image.open(rp)
    row["ours_wh"] = "%dx%d" % oi.size
    row["ref_wh"] = "%dx%d" % ri.size
    row["aspect_ours"] = round(oi.size[0] / oi.size[1], 3)
    row["aspect_ref"] = round(ri.size[0] / ri.size[1], 3)
    row["aspect_ratio"] = round(min(row["aspect_ours"], row["aspect_ref"]) /
                                max(row["aspect_ours"], row["aspect_ref"]), 3)
    sa, sb = structural(op), structural(rp)
    if sa is None or sb is None:
        row["structural_r"] = "na(blank)"
    else:
        row["structural_r"] = round(float(np.dot(sa.ravel(), sb.ravel())), 3)
    ma, mb = ink_mask(op), ink_mask(rp)
    inter = np.logical_and(ma, mb).sum()
    union = np.logical_or(ma, mb).sum()
    blank_ref = mb.mean() < 0.01
    row["ink_iou"] = "na(blank)" if blank_ref else (round(float(inter / union), 3) if union else 0.0)
    row["ink_ours"] = round(float(ma.mean()), 3)
    row["ink_ref"] = round(float(mb.mean()), 3)
    no, co = colours(op)
    nr, cr = colours(rp)
    row["colours_ours"] = no
    row["colours_ref"] = nr
    row["top_colours_ours"] = ";".join("%s %.3f" % (k, v) for k, v in list(co.items())[:6])
    row["top_colours_ref"] = ";".join("%s %.3f" % (k, v) for k, v in list(cr.items())[:6])
    return row


def main():
    run_root, ref_root, out_tsv = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    sub = "figure"
    if "--figs" in sys.argv:
        sub = sys.argv[sys.argv.index("--figs") + 1]
    rows = []
    for ours, ref26, refhtml, cell, what, kind in PAIRS:
        op = run_root / sub / ours
        rp26 = ref_root / "teacher2026" / ref26
        rph = ref_root / "teacherhtml" / refhtml if refhtml else None
        row = {
            "figure": ours,
            "tutorial_cell": cell,
            "tutorial_call": what,
            "kind": kind,
            "ours_exists": op.exists(),
            "ref2026_file": ref26,
            "ref2026_exists": rp26.exists(),
            "refhtml_file": refhtml or "",
        }
        if op.exists() and rp26.exists():
            for k, v in metrics(op, rp26).items():
                row["vs2026_" + k] = v
        if op.exists() and rph is not None and rph.exists():
            for k, v in metrics(op, rph).items():
                row["vshtml_" + k] = v
        rows.append(row)

    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("%-46s %-5s %-5s %-6s %-6s %-6s %-6s" %
          ("figure", "2026r", "2026iou", "herr", "hiou", "ourA", "refA"))
    for r in rows:
        print("%-46s %-5s %-7s %-6s %-6s %-6s %-6s" % (
            r["figure"],
            r.get("vs2026_structural_r", "-"),
            r.get("vs2026_ink_iou", "-"),
            r.get("vshtml_structural_r", "-"),
            r.get("vshtml_ink_iou", "-"),
            r.get("vs2026_aspect_ours", "-"),
            r.get("vs2026_aspect_ref", "-")))
    print("\n-> %s" % out_tsv)


if __name__ == "__main__":
    main()
