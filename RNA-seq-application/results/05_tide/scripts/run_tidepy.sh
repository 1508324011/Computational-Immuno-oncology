#!/usr/bin/env bash
## =============================================================================
## Lesson 3 — RNA-seq application: Part 4 TIDE — shell step (code block 27)
## Source: RNA-seq-application.html
## Deviations (pre-registered, additive only):
##   - run via sbatch-free login-node serial (RNA-seq discipline #52, not WES sbatch)
##   - output path made concrete (tutorial placeholder: /path_to_your_folder/)
##     and written into this step's results dir
## =============================================================================
set -eu
source /lustre1/share/miniconda3/etc/profile.d/conda.sh
conda activate rnaimmune

WORKDIR="/gpfs1/home/bjx131_pkuhpc/class/RNA-seq-application/results/05_tide"
EXPR="/lustre1/share/rnaimmune/data/expression/VanAllen2015_CTLA4_Melanoma_RNASeq"
OUT="${WORKDIR}/res_tide.txt"

cd "$WORKDIR"
echo "[tidepy] start $(date -Is)"
echo "[tidepy] cmd: tidepy ${EXPR} -o ${OUT} -c Melanoma --force_normalize"
# code block 27 (verbatim except /path_to_your_folder/ -> concrete workdir path)
tidepy "${EXPR}" -o "${OUT}" -c Melanoma --force_normalize 2>&1
echo "[tidepy] done exit=$? $(date -Is)"
echo "[tidepy] output:"
ls -la "${OUT}"
echo "[tidepy] head:"
head -3 "${OUT}"
echo "[tidepy] rows:"
wc -l "${OUT}"
