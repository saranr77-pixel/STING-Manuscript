"""Heatmap of STING1 CpG methylation across TCGA-OV tumors.

Input: data/TCGA_OV_STING1_CpG_methylation_xena.xlsx (Xena HM450 export).
Only the samples carrying HM450 data are plotted; beta 0 = unmethylated,
1 = fully methylated. Samples are ordered by mean beta across the probes.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

REGION = {}   # probe -> region label (TSS200 etc.); fill from the manifest

df = pd.read_excel("data/TCGA_OV_STING1_CpG_methylation_xena.xlsx")
probes = [c for c in df.columns if c.startswith("cg")]
df[probes] = df[probes].apply(pd.to_numeric, errors="coerce")
df = df[df[probes].notna().any(axis=1)].copy()
df["mean_beta"] = df[probes].mean(axis=1)
df = df.sort_values("mean_beta", ascending=False)
mat = df[probes].to_numpy().T           # probes x samples

cmap = LinearSegmentedColormap.from_list("meth", ["#F7F7F7", "#F2C0B4", "#C0392B"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, ax = plt.subplots(figsize=(7.4, 3.9))
im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)

ax.set_yticks(range(len(probes)))
ax.set_yticklabels([f"{p}  {REGION[p]}" if REGION.get(p) else p for p in probes],
                   fontsize=9)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df["sample"], rotation=90, fontsize=6.5)
ax.set_xlabel(f"TCGA-OV primary tumors (n = {len(df)}), ordered by mean methylation")
ax.set_title("STING1 (TMEM173) CpG methylation", fontsize=11, pad=8)
for x in np.arange(0.5, len(df)):
    ax.axvline(x, color="white", linewidth=0.6)
for y in np.arange(0.5, len(probes)):
    ax.axhline(y, color="white", linewidth=0.6)
cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.035)
cb.set_label("$\\beta$ value", fontsize=9)
cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/Fig_STING1_CpG_heatmap.{ext}", dpi=300)
print(f"{len(probes)} probes x {len(df)} samples")
