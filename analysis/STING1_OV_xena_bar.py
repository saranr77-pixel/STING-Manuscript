"""STING1 (TMEM173) in TCGA-OV tumors vs GTEx normal ovary -- bar graph.

Input: data/Normal_vs_tumor_UCSC_xena.xlsx, exported from the UCSC Xena browser
(TCGA TARGET GTEx cohort). Values are Xena's RSEM norm_count, log2(x+1).
Bars = mean +/- SD, every sample shown; two-sided Mann-Whitney U test.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

rng = np.random.default_rng(1)
df = pd.read_excel("data/Normal_vs_tumor_UCSC_xena.xlsx", header=2)
df["expr"] = pd.to_numeric(df["TMEM173"])

groups = [("Normal ovary\n(GTEx)", "Normal Tissue", "#2E6FBF"),
          ("Ovarian cancer\n(TCGA-OV)", "Primary Tumor", "#D6453D")]
vals = [df.loc[df["_sample_type"] == st, "expr"].to_numpy() for _, st, _ in groups]

u = stats.mannwhitneyu(vals[0], vals[1], alternative="two-sided")
excl = vals[1][vals[1] > 0]       # sensitivity: drop the one zero-count tumor
u_excl = stats.mannwhitneyu(vals[0], excl, alternative="two-sided")
diff = vals[1].mean() - vals[0].mean()

rows = []
for (name, _, _), v in zip(groups, vals):
    rows.append(dict(group=name.replace("\n", " "), n=len(v), mean=v.mean(),
                     sd=v.std(ddof=1), median=np.median(v),
                     q1=np.percentile(v, 25), q3=np.percentile(v, 75)))
summ = pd.DataFrame(rows)
summ["log2FC_tumor_vs_normal"] = diff
summ["MannWhitney_p"] = u.pvalue
summ["MannWhitney_p_excl_zero"] = u_excl.pvalue
summ.to_csv("figures/STING1_OV_xena_stats.csv", index=False)
print(summ.to_string(index=False))

def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < .01 else "*" if p < .05 else "ns"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig, ax = plt.subplots(figsize=(4.0, 4.6))
for i, ((name, _, col), v) in enumerate(zip(groups, vals)):
    ax.bar(i, v.mean(), width=0.62, color=col, alpha=0.85, edgecolor="black",
           linewidth=0.8, zorder=1)
    ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=9, facecolor=col,
               edgecolor="black", linewidth=0.3, alpha=0.55, zorder=2)
    ax.errorbar(i, v.mean(), yerr=v.std(ddof=1), color="black", capsize=6,
                linewidth=1.2, capthick=1.2, zorder=3)

top = max(v.max() for v in vals)
y = top + 0.5
ax.plot([0, 0, 1, 1], [y - 0.2, y, y, y - 0.2], color="black", linewidth=1)
ax.text(0.5, y + 0.1, f"{stars(u.pvalue)}\np = {u.pvalue:.3f}", ha="center",
        va="bottom", fontsize=10)

ax.set_xticks([0, 1])
ax.set_xticklabels([f"{n}\nn = {len(v)}" for (n, _, _), v in zip(groups, vals)])
ax.set_ylabel("STING1 (TMEM173) expression\nlog$_2$(RSEM norm_count + 1)")
ax.set_ylim(0, y + 1.6)
ax.set_xlim(-0.6, 1.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/Fig_STING1_OV_vs_normal_bar.{ext}", dpi=300)
