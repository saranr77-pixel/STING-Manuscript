"""STING1 (TMEM173) tumor (TCGA) vs normal (GTEx) across tissues -- bar graph.

Input: data/Pancancer_normal_vs_tumor_UCSC_xena.xlsx, one sheet per primary
site, exported from the UCSC Xena browser (TCGA TARGET GTEx cohort).
Values are Xena's RSEM norm_count, log2(x+1).
Bars = mean +/- SD, every sample shown. Two-sided Mann-Whitney U per tissue,
Benjamini-Hochberg FDR across tissues. Samples with no value are dropped;
samples with a value of exactly 0 (failed libraries) are dropped and counted.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.patches import Patch

SRC = "data/Pancancer_normal_vs_tumor_UCSC_xena.xlsx"
NORMAL, TUMOR = "#2E6FBF", "#D6453D"
rng = np.random.default_rng(1)

frames = []
for sheet in pd.ExcelFile(SRC).sheet_names:
    raw = pd.read_excel(SRC, sheet_name=sheet, header=None)
    h = raw.index[raw[0] == "sample"][0]
    d = raw.iloc[h + 1:].copy()
    d.columns = raw.iloc[h].tolist()
    frames.append(d)
df = pd.concat(frames, ignore_index=True)
df["expr"] = pd.to_numeric(df["TMEM173"], errors="coerce")
df["tissue"] = df["_primary_site"]
df["group"] = df["_sample_type"].map({"Normal Tissue": "Normal", "Primary Tumor": "Tumor"})

n_missing = df.groupby(["tissue", "group"])["expr"].apply(lambda s: s.isna().sum())
n_zero = df.groupby(["tissue", "group"])["expr"].apply(lambda s: (s == 0).sum())
df = df[df["expr"].notna() & (df["expr"] > 0)]

rows = []
for t, g in df.groupby("tissue"):
    nv, tv = g.loc[g.group == "Normal", "expr"], g.loc[g.group == "Tumor", "expr"]
    rows.append(dict(tissue=t, n_normal=len(nv), n_tumor=len(tv),
                     mean_normal=nv.mean(), sd_normal=nv.std(),
                     mean_tumor=tv.mean(), sd_tumor=tv.std(),
                     median_normal=nv.median(), median_tumor=tv.median(),
                     diff_log2=tv.mean() - nv.mean(),
                     p=stats.mannwhitneyu(tv, nv, alternative="two-sided").pvalue,
                     zeros_removed_normal=n_zero[(t, "Normal")],
                     zeros_removed_tumor=n_zero[(t, "Tumor")],
                     missing_normal=n_missing[(t, "Normal")]))
res = pd.DataFrame(rows)
res["FDR"] = stats.false_discovery_control(res["p"], method="bh")
res["direction"] = np.where(res["FDR"] >= .05, "ns",
                            np.where(res["diff_log2"] < 0, "lower in tumor", "higher in tumor"))
res = res.sort_values("tissue").reset_index(drop=True)
res.to_csv("figures/STING1_pancancer_xena_stats.csv", index=False)
print(res.round(4).to_string(index=False))

def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < .01 else "*" if p < .05 else "ns"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, ax = plt.subplots(figsize=(11, 4.8))
w = 0.38
top = df["expr"].max()
for i, r in res.iterrows():
    g = df[df.tissue == r.tissue]
    for off, grp, col in ((-w / 2, "Normal", NORMAL), (w / 2, "Tumor", TUMOR)):
        v = g.loc[g.group == grp, "expr"].to_numpy()
        ax.bar(i + off, v.mean(), width=w * 0.92, color=col, alpha=0.85,
               edgecolor="black", linewidth=0.7, zorder=1,
               )
        ax.scatter(i + off + rng.uniform(-w * 0.3, w * 0.3, len(v)), v, s=3,
                   facecolor=col, edgecolor="black", linewidth=0.15, alpha=0.35,
                   zorder=2)
        ax.errorbar(i + off, v.mean(), yerr=v.std(ddof=1), color="black",
                    capsize=3, linewidth=1, capthick=1, zorder=3)
    y = top + 0.35
    ax.plot([i - w / 2, i - w / 2, i + w / 2, i + w / 2], [y - .15, y, y, y - .15],
            color="black", linewidth=0.8)
    ax.text(i, y + 0.05, stars(r.FDR), ha="center", va="bottom", fontsize=9)

ax.set_xticks(range(len(res)))
ax.set_xticklabels([f"{t}\nN={a}\nT={b}" for t, a, b in
                    zip(res.tissue, res.n_normal, res.n_tumor)])
for lab in ax.get_xticklabels():
    if lab.get_text().startswith("Ovary"):
        lab.set_fontweight("bold")
ax.set_ylabel("STING1 (TMEM173) expression\nlog$_2$(RSEM norm_count + 1)")
ax.set_ylim(0, top + 1.3)
ax.set_xlim(-0.6, len(res) - 0.4)
ax.legend(handles=[Patch(facecolor=NORMAL, edgecolor="black", label="Normal (GTEx)"),
                   Patch(facecolor=TUMOR, edgecolor="black", label="Tumor (TCGA)")],
          frameon=False, ncol=2,
          loc="upper right", bbox_to_anchor=(1, 1.1))
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/Fig_STING1_pancancer_bar.{ext}", dpi=300)
