"""STING1 (TMEM173) CpG methylation by probe, ordered along the gene.

Layout follows Kitajima/Xia-style promoter methylation panels: one horizontal
box per CpG probe, probes ordered by genomic position (STING1 is on the minus
strand of chr5, so the TSS is at the HIGH-coordinate end, plotted at the top),
with a gene-region annotation bar on the left.

Data: data/TCGA_OV_STING1_CpG_methylation_xena.xlsx (Xena HM450 export).
LIMITATION: only 7 of the 11 STING1 probes are in this export, only 10 TCGA-OV
tumors carry HM450 data, and the export has no usable normal samples. Swap in a
probe-level dataset with tumours AND normals to produce the publication panel.

Coordinates and gene-region assignments transcribed from the reference figure --
VERIFY each against the Illumina manifest before publishing.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# probe -> (chr5 coordinate, gene region). Ordered 5'->3' for STING1 (minus strand).
ANNOT = [
    ("cg21618694", 138863439, "TSS1500"),
    ("cg10813908", 138863041, "TSS1500"),
    ("cg00070715", 138862796, "TSS1500"),
    ("cg16983159", 138862441, "TSS200"),
    ("cg23255964", 138862353, "TSS200"),
    ("cg16532438", 138861941, "5'UTR"),
    ("cg14655316", 138861855, "1stExon"),
    ("cg04560810", 138861847, "1stExon"),
    ("cg03317505", 138861832, "1stExon"),
    ("cg04232128", 138861241, "Body"),
    ("cg01938023", 138855699, "3'UTR"),
]
REGION_COLOR = {"TSS1500": "#E03B3B", "TSS200": "#2E86C1", "5'UTR": "#27AE60",
                "1stExon": "#E67E22", "Body": "#C9B037", "3'UTR": "#7B4A12"}
TUMOR = "#1F3F8F"

df = pd.read_excel("data/TCGA_OV_STING1_CpG_methylation_xena.xlsx")
have = [c for c in df.columns if c.startswith("cg")]
df[have] = df[have].apply(pd.to_numeric, errors="coerce")
df = df[df[have].notna().any(axis=1)]

annot = [a for a in ANNOT if a[0] in have]
missing = [a[0] for a in ANNOT if a[0] not in have]
print(f"plotted {len(annot)} probes, n={len(df)} tumors; no data for: {missing}")

fig, (cax, ax) = plt.subplots(
    1, 2, figsize=(8.6, 4.6), gridspec_kw={"width_ratios": [0.05, 1], "wspace": 0.02})
plt.rcParams.update({"font.family": "DejaVu Sans"})

ys = np.arange(len(annot))[::-1]          # first probe (5') at the top
for y, (probe, coord, region) in zip(ys, annot):
    v = df[probe].dropna().to_numpy()
    ax.boxplot(v, positions=[y], vert=False, widths=0.55, whis=1.5,
               patch_artist=True, showfliers=True,
               boxprops=dict(facecolor=TUMOR, color=TUMOR),
               medianprops=dict(color="white", linewidth=1.2),
               whiskerprops=dict(color=TUMOR, linewidth=1),
               capprops=dict(color=TUMOR, linewidth=1),
               flierprops=dict(marker="o", markersize=3, markerfacecolor="none",
                               markeredgecolor=TUMOR, markeredgewidth=0.6))
    cax.add_patch(plt.Rectangle((0, y - 0.5), 1, 1,
                                facecolor=REGION_COLOR[region], edgecolor="none"))

cax.set_xlim(0, 1); cax.set_ylim(-0.5, len(annot) - 0.5)
cax.set_xticks([]); cax.set_yticks(ys)
cax.set_yticklabels([f"{c:,}" for _, c, _ in annot], fontsize=7.5)
cax.set_xlabel("Gene", fontsize=8)
for s in cax.spines.values():
    s.set_visible(False)

ax.set_ylim(-0.5, len(annot) - 0.5)
ax.set_yticks(ys)
ax.set_yticklabels([p for p, _, _ in annot], fontsize=8.5)
ax.yaxis.tick_right()
ax.set_xlim(0, 1)
ax.set_xlabel("$\\beta$ value")
ax.set_title(f"STING1 (TMEM173) — TCGA-OV primary tumors (n = {len(df)})",
             fontsize=10, pad=8)
ax.grid(axis="x", color="0.9", linewidth=0.6)
ax.set_axisbelow(True)
for s in ("top", "left", "right"):
    ax.spines[s].set_visible(False)

order = ["TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR"]
present = [r for r in order if any(a[2] == r for a in annot)]
fig.legend(handles=[Patch(facecolor=REGION_COLOR[r], label=r) for r in present],
           loc="lower center", ncol=len(present), frameon=False, fontsize=8,
           bbox_to_anchor=(0.55, -0.09))
fig.tight_layout(rect=[0, 0.10, 1, 1])
for ext in ("pdf", "png"):
    fig.savefig(f"figures/Fig_STING1_CpG_bygene.{ext}", dpi=300, bbox_inches="tight")
