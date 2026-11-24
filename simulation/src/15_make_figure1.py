"""Figure 1: sensitivity and false-flag proportion of the primary rule by k and expected-count stratum (pooled over replicates)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
O = HERE.parent / "outputs"
d = pd.read_csv(O / "grid_summary.csv")
d = d[(d.rule == "post_tau1.10_c0.95") & (d.k > 0)]
strata = [("<20", "Fewer than 20", "#9ecae1", "o", ":"), ("20-49", "20-49", "#6baed6", "s", "-."),
          ("50-199", "50-199", "#3182bd", "^", "--"), (">=200", "200 or more", "#08306b", "D", "-")]
ks = sorted(d.k.unique())
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=300, sharex=True)
for a, col, title, ymax in [(ax[0], "sens_pooled", "A. Sensitivity", 0.8), (ax[1], "false_flag_pooled", "B. Proportion of flags that are false", 0.04)]:
    for lab, name, c, mk, ls in strata:
        x = d[d.stratum == lab].sort_values("k")
        a.plot(x.k, x[col], color=c, marker=mk, linestyle=ls, linewidth=1.6, markersize=5, markeredgecolor="#08306b", markeredgewidth=0.5, label=name)
    a.set_xlabel("Scale of generating heterogeneity, k")
    a.set_xticks(ks)
    a.set_ylim(0, ymax)
    a.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0 if col == "sens_pooled" else 1))
    a.set_title(title, loc="left", fontsize=9, fontweight="bold")
    a.grid(axis="y", color="#dddddd", linewidth=0.6)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
ax[0].legend(title="Expected deaths", frameon=False, loc="upper left", fontsize=7.5, title_fontsize=7.5)
fig.tight_layout()
fig.savefig(HERE.parent / "outputs" / "Figure_1.png", dpi=300)
fig.savefig(HERE.parent / "outputs" / "Figure_1.pdf")
m = d.groupby("stratum")[["sens_rep_mcse", "false_flag_rep_mcse"]].max()
print(m.round(4))
print("max false-flag value:", d.false_flag_pooled.max().round(4))
