"""STING1 promoter methylation (cg16983159) vs STING1 mRNA in TCGA-OV.

Input: data/cbioportal_STING1_meth_vs_expr.xlsx, exported from cBioPortal.
  x = TMEM173 cg16983159 beta value (HM27 and HM450 merge)
  y = STING1 mRNA expression z-score (log RNA Seq V2 RSEM), relative to all samples

Reports Pearson and Spearman correlations with a fitted line and 95% CI band.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

d = pd.read_excel("data/cbioportal_STING1_meth_vs_expr.xlsx", header=1)
d.columns = ["sample", "meth", "expr"]
d[["meth", "expr"]] = d[["meth", "expr"]].apply(pd.to_numeric, errors="coerce")
d = d.dropna(subset=["meth", "expr"])
n = len(d)

r, rp = stats.pearsonr(d.meth, d.expr)
rho, sp = stats.spearmanr(d.meth, d.expr)
fit = stats.linregress(d.meth, d.expr)
pd.DataFrame([dict(n=n, pearson_r=r, pearson_p=rp, r_squared=r**2,
                   spearman_rho=rho, spearman_p=sp,
                   slope=fit.slope, intercept=fit.intercept,
                   n_beta_gt_0_5=int((d.meth > 0.5).sum()))]
             ).to_csv("figures/STING1_meth_vs_expr_stats.csv", index=False)
print(f"n={n}  Pearson r={r:.3f} p={rp:.2e} R2={r**2:.3f}  "
      f"Spearman rho={rho:.3f} p={sp:.2e}")

# 95% CI band for the fitted line
xs = np.linspace(d.meth.min(), d.meth.max(), 200)
ys = fit.intercept + fit.slope * xs
resid = d.expr - (fit.intercept + fit.slope * d.meth)
se = np.sqrt((resid**2).sum() / (n - 2)) * np.sqrt(
    1 / n + (xs - d.meth.mean())**2 / ((d.meth - d.meth.mean())**2).sum())
tcrit = stats.t.ppf(0.975, n - 2)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig, ax = plt.subplots(figsize=(5.2, 4.8))
ax.scatter(d.meth, d.expr, s=18, facecolor="#2E6FBF", edgecolor="black",
           linewidth=0.3, alpha=0.55, zorder=2)
ax.fill_between(xs, ys - tcrit * se, ys + tcrit * se, color="#C0392B",
                alpha=0.18, zorder=3, linewidth=0)
ax.plot(xs, ys, color="#C0392B", linewidth=2, zorder=4)

txt = (f"$r$ = {r:.2f}  ($p$ = {rp:.1e})\n"
       f"$\\rho$ = {rho:.2f}  ($p$ = {sp:.1e})\n"
       f"$R^2$ = {r**2:.2f}   $n$ = {n}")
ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top",
        fontsize=9.5, bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                                edgecolor="lightgrey"))
ax.set_xlabel("STING1 (TMEM173) methylation, cg16983159\n"
              "$\\beta$ value (HM27/HM450 merge)")
ax.set_ylabel("STING1 mRNA expression\n($z$-score, log RNA Seq V2 RSEM)")
ax.set_xlim(0, 1)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/Fig_STING1_meth_vs_expr.{ext}", dpi=300)
