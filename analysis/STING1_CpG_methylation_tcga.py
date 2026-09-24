"""STING1 (TMEM173) CpG methylation across TCGA-OV tumors.

Input: data/TCGA_OV_STING1_CpG_methylation_xena.xlsx, exported from the UCSC
Xena browser. Beta values (0 = unmethylated, 1 = fully methylated).

NOTE: this export is the HM450 dataset, in which only 10 TCGA-OV primary
tumors carry data; the remaining rows (incl. all Solid Tissue Normal) are
empty. Re-export the HM27/HM450 merged dataset for the full cohort.

REGION labels below are placeholders keyed by probe ID -- fill them in from
the Illumina manifest (UCSC_RefGene_Group) before publishing.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "data/TCGA_OV_STING1_CpG_methylation_xena.xlsx"
GENE = "STING1 (TMEM173)"
# probe -> region relative to the STING1 TSS. Verify against the manifest.
REGION = {
    "cg16983159": "",
    "cg23255964": "",
    "cg16532438": "",
    "cg04560810": "",
    "cg03317505": "",
    "cg04232128": "",
    "cg01938023": "",
}
HIGHLIGHT = {"cg16983159"}          # probe(s) to emphasise
METH, UNMETH = "#C0392B", "#2E6FBF"  # >0.5 methylated / <=0.5 unmethylated
rng = np.random.default_rng(0)

df = pd.read_excel(SRC)
probes = [c for c in df.columns if c.startswith("cg")]
df[probes] = df[probes].apply(pd.to_numeric, errors="coerce")
df = df[df[probes].notna().any(axis=1)]
print(df["sample_type"].value_counts().to_string())

summ = pd.DataFrame({
    "probe": probes,
    "region": [REGION.get(p, "") for p in probes],
    "n": [df[p].notna().sum() for p in probes],
    "mean_beta": [df[p].mean() for p in probes],
    "sd": [df[p].std(ddof=1) for p in probes],
    "median_beta": [df[p].median() for p in probes],
    "min": [df[p].min() for p in probes],
    "max": [df[p].max() for p in probes],
    "pct_samples_beta_gt_0.5": [(df[p] > 0.5).mean() * 100 for p in probes],
})
summ.to_csv("figures/STING1_CpG_methylation_tcga_stats.csv", index=False)
print(summ.round(3).to_string(index=False))

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig, ax = plt.subplots(figsize=(7.6, 4.8))
ax.axhspan(0.5, 1.0, color="#C0392B", alpha=0.05, zorder=0)
ax.axhline(0.5, color="grey", linestyle="--", linewidth=1, zorder=1)
ax.text(len(probes) - 0.45, 0.515, "methylated ($\\beta$ > 0.5)", ha="right",
        va="bottom", fontsize=9, color="grey")

for i, p in enumerate(probes):
    v = df[p].dropna().to_numpy()
    col = METH if v.mean() > 0.5 else UNMETH
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=26, facecolor=col,
               edgecolor="black", linewidth=0.4, alpha=0.75, zorder=3)
    ax.plot([i - 0.28, i + 0.28], [v.mean()] * 2, color="black",
            linewidth=1.8, zorder=4)
    ax.errorbar(i, v.mean(), yerr=v.std(ddof=1), color="black", capsize=4,
                linewidth=1, zorder=4)

labels = []
for p in probes:
    lab = f"$\\bf{{{p}}}$" if p in HIGHLIGHT else p
    labels.append(f"{lab}\n{REGION[p]}" if REGION.get(p) else lab)
ax.set_xticks(range(len(probes)))
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("DNA methylation ($\\beta$ value)")
ax.set_title(f"{GENE} CpG methylation in TCGA-OV primary tumors "
             f"(n = {len(df)})", fontsize=11)
ax.set_ylim(0, 1.0)
ax.set_xlim(-0.6, len(probes) - 0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/Fig_STING1_CpG_methylation_tcga.{ext}", dpi=300)
