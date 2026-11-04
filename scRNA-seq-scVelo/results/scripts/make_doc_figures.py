#!/usr/bin/env python3
"""Build the composite / derived figures used by the two corrected docs and the results guide.

Three kinds of figure are produced here, all into <module>/ref/figures/:

  1. comparison panels   teacher reference on top, our as-written render underneath, so the
                         per-panel correlation quoted in the text can be eyeballed.
  2. derived bar charts  the figure-fidelity table and the parameter-sensitivity table turned
                         into pictures (numbers identical to results/07_figure_comparison/*.tsv
                         and results/06_corrected/sensitivity_*.csv).
  3. pass-throughs      the corrected-pass panels the prose refers to, copied next to the docs
                         so the markdown/HTML/PDF stay self-contained.

Usage:  python make_doc_figures.py <run_root> <module_root>
"""
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

RUN = Path(sys.argv[1])
MOD = Path(sys.argv[2])
FIG = MOD / "ref" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------- 1. comparison panels
def stack(ref_png, ours_png, out, caption_ref, caption_ours, gap=16, bg=(255, 255, 255)):
    """Stack two PNGs vertically at a common width; labels are drawn in the caption, not here."""
    a, b = Image.open(ref_png).convert("RGB"), Image.open(ours_png).convert("RGB")
    w = max(a.width, b.width)
    a = a.resize((w, round(a.height * w / a.width)), Image.LANCZOS)
    b = b.resize((w, round(b.height * w / b.width)), Image.LANCZOS)
    canvas = Image.new("RGB", (w, a.height + gap + b.height), bg)
    canvas.paste(a, (0, 0))
    canvas.paste(b, (0, a.height + gap))
    canvas.save(out, optimize=True)
    print("wrote", out.name, canvas.size)


stack(MOD / "02_velocity_embedding_stream.png", MOD / "results/02_velocity/02_velocity_embedding_stream.png",
      FIG / "compare_stream.png", "reference", "as-written")

stack(MOD / "08_cellcycle_gradients.png", MOD / "results/05_cellcycle/08_cellcycle_gradients.png",
      FIG / "compare_cellcycle.png", "reference", "as-written")

stack(MOD / "04_phaseportrait_markers.png",
      MOD / "results/06_corrected/corrected_05_phaseportrait_markers_dynamical.png",
      FIG / "compare_phase_dynamical.png", "reference", "dynamical fits")

# ----------------------------------------------------------------- 2. derived bar charts
tsv = pd.read_csv(MOD / "results/07_figure_comparison/figure_comparison_aswritten.tsv", sep="\t")
t = tsv.dropna(subset=["struct_corr"]).copy()
t["short"] = t["figure"].str.replace(r"\.png$", "", regex=True)
fig, ax = plt.subplots(figsize=(7.2, 3.4), dpi=300)
y = np.arange(len(t))[::-1]
ax.barh(y, t["struct_corr"], color=["#c0392b" if c < 0.9 else "#2e6da4" for c in t["struct_corr"]])
ax.set_yticks(y)
ax.set_yticklabels(t["short"], fontsize=6.5)
ax.set_xlim(0, 1.05)
ax.set_xlabel("structure correlation against the reference PNG", fontsize=8)
ax.axvline(0.9, color="grey", lw=0.6, ls=":")
for yy, c, i in zip(y, t["struct_corr"], t["ink_iou"]):
    ax.text(c + 0.005, yy, f"{c:.3f} / IoU {i:.2f}", va="center", fontsize=6)
ax.set_title("as-written run vs the 11 figures shipped with the course", fontsize=9)
ax.tick_params(axis="x", labelsize=7)
fig.tight_layout()
fig.savefig(FIG / "figure_fidelity_aswritten.png")
plt.close(fig)

sens = pd.read_csv(MOD / "results/06_corrected/sensitivity_graph_and_hvg.csv")
fig, axes = plt.subplots(1, 2, figsize=(7.6, 2.9), dpi=300)
lab = sens["variant"].tolist()
jac = "top20_jaccard_vs_stochastic_main"
axes[0].barh(range(len(sens))[::-1], sens[jac], color="#2e6da4")
axes[0].set_yticks(range(len(sens))[::-1])
axes[0].set_yticklabels(lab, fontsize=6.5)
axes[0].set_xlim(0, 1.05)
axes[0].set_xlabel("top-20 velocity genes: Jaccard vs the default run", fontsize=7)
axes[0].axvline(1.0, color="grey", lw=0.6, ls=":")
axes[1].barh(range(len(sens))[::-1], sens["n_velocity_genes"], color="#8e44ad")
axes[1].set_yticks(range(len(sens))[::-1])
axes[1].set_yticklabels([])
axes[1].set_xlabel("genes called velocity genes", fontsize=7)
for i, (v, n) in enumerate(zip(sens[jac], sens["n_velocity_genes"])):
    axes[0].text(min(v + 0.01, 0.96), len(sens) - 1 - i, f"{v:.2f}", va="center", fontsize=6)
    axes[1].text(n + 8, len(sens) - 1 - i, f"{n}", va="center", fontsize=6)
axes[1].set_xlim(0, max(sens["n_velocity_genes"]) * 1.25)
fig.suptitle("what the tutorial's unexplained defaults cost (dynamical model, stochastic velocity)",
             fontsize=8.5)
fig.tight_layout()
fig.savefig(FIG / "sensitivity_defaults.png")
plt.close(fig)

# ----------------------------------------------------------------- 3. pass-throughs
for p in ["corrected_01_proportions.png", "corrected_02_velocity_embedding_stream.png",
          "corrected_04_phaseportrait_markers_stochastic.png", "corrected_05_phaseportrait_markers_dynamical.png",
          "corrected_09_cellcycle_gradients_original_scores.png",
          "corrected_11_phaseportrait_hells_top2a_dynamical.png",
          "corrected_12_latent_time_and_confidence.png",
          "corrected_13_velocity_stream_latent_time.png"]:
    shutil.copy2(MOD / "results/06_corrected" / p, FIG / p)
print("copied", 8, "corrected panels")

# theory figures the notebook embeds, archived under ref/figures/ with a theory_ prefix
for name in ["stability", "func", "intron-exon", "sss", "ssd", "sst", "ss1", "scvelo", "label"]:
    src = RUN / "theory_img" / f"{name}.png"
    if src.exists():
        shutil.copy2(src, FIG / f"theory_{name}.png")
print("done")
