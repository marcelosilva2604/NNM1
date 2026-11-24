"""Figure 3 (operating points of the rule family and comparators) and Figure 1 (map of the registry application).

Figure 2 reads outputs/grid_summary.csv. Figure 3 needs the municipal mesh: set MUNI_GEOJSON and STATES_GEOJSON to the polygon files
(feature id = 7-digit IBGE code, first 6 digits = analysis code) and PRIMARY_MUNI to model/muni_bym2_rsr_full.csv.
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
O = HERE.parent / "outputs"
SUB = O
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})

# ---------------- Figure 3 (operating points)
d = pd.read_csv(O / "grid_summary.csv")
d = d[d.stratum == "ALL"]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), dpi=300)
fam = [("post_tau1.10_c0.99", "c = 0.99"), ("post_tau1.10_c0.95", "c = 0.95"), ("post_tau1.10_c0.90", "c = 0.90"), ("post_tau1.10_c0.80", "c = 0.80")]
other_tau = ["post_tau1.05_c0.95", "post_tau1.20_c0.95"]
comp = [("C1_median_gt1.10", "C1"), ("C2_poisson_p05", "C2"), ("C3_poisson_BH05", "C3")]
cfg = {1.0: dict(ylim=(-3, 50), xlim=(-1, 34), lab="A. k = 1.0", off={"post_tau1.10_c0.99": (8, 0), "post_tau1.10_c0.95": (8, -4), "post_tau1.10_c0.90": (8, 0), "post_tau1.10_c0.80": (8, 2), "C1_median_gt1.10": (7, 3), "C2_poisson_p05": (7, 3), "C3_poisson_BH05": (7, -8)}),
       2.0: dict(ylim=(-3, 75), xlim=(-1, 26), lab="B. k = 2.0", off={"post_tau1.10_c0.99": (8, 0), "post_tau1.10_c0.95": (8, -5), "post_tau1.10_c0.90": (8, 6), "post_tau1.10_c0.80": (8, 2), "C1_median_gt1.10": (7, 3), "C2_poisson_p05": (7, -7), "C3_poisson_BH05": (7, -8)})}
for ax, k in [(axs[0], 1.0), (axs[1], 2.0)]:
    c = cfg[k]
    x = d[d.k == k].set_index("rule")
    fx = lambda r: 100 * x.loc[r, "false_flag_pooled"]
    fy = lambda r: 100 * x.loc[r, "sens_pooled"]
    ax.plot([fx(r) for r, _ in fam], [fy(r) for r, _ in fam], color="#08306b", marker="o", markersize=5, linewidth=1.4, label="P(RR > 1.10) ≥ c", zorder=3)
    for r, t in fam:
        ax.annotate(t, (fx(r), fy(r)), textcoords="offset points", xytext=c["off"][r], fontsize=6.3, color="#08306b", va="center")
    ax.plot([fx(r) for r in other_tau], [fy(r) for r in other_tau], marker="^", markersize=5, markerfacecolor="white", markeredgecolor="#08306b", linestyle="none", label="τ = 1.05 and 1.20, c = 0.95", zorder=3)
    ax.plot([fx(r) for r, _ in comp], [fy(r) for r, _ in comp], marker="s", markersize=5, color="#7f7f7f", linestyle="none", label="Comparators C1, C2, C3", zorder=3)
    for r, t in comp:
        ax.annotate(t, (fx(r), fy(r)), textcoords="offset points", xytext=c["off"][r], fontsize=6.8, color="#404040", va="center")
    ax.set_xlabel("False flags, % of flagged")
    ax.set_ylabel("Sensitivity, %")
    ax.set_xlim(*c["xlim"]); ax.set_ylim(*c["ylim"])
    ax.set_title(c["lab"], loc="left", fontsize=9, fontweight="bold")
    ax.grid(color="#dddddd", linewidth=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
axs[0].legend(frameon=False, loc="center right", fontsize=6.8)
fig.tight_layout()
fig.savefig(SUB / "Figure_3.png", dpi=300)
fig.savefig(SUB / "Figure_3.pdf")
plt.close(fig)
print("figure 3 (operating points) ok")

# ---------------- Figure 1 (map)
gj = os.environ.get("MUNI_GEOJSON")
if gj:
    import geopandas as gpd
    from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
    from matplotlib.patches import Patch
    CRS = 5880
    m = pd.read_csv(os.environ["PRIMARY_MUNI"], dtype={"CODMUNRES": str})
    munis = gpd.read_file(gj).to_crs(CRS)
    munis["CODMUNRES"] = munis["id"].astype(str).str[:6]
    g = munis.merge(m, on="CODMUNRES", how="inner")
    states = gpd.read_file(os.environ["STATES_GEOJSON"]).to_crs(CRS)

    def scale_and_north(ax, km=1000):
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); W, H = x1 - x0, y1 - y0
        L = km * 1000.0
        bx, by = x1 - 0.04 * W - L, y0 + 0.05 * H
        ax.plot([bx, bx + L], [by, by], color="black", lw=2.0, solid_capstyle="butt", zorder=10)
        ax.text(bx + L / 2, by + 0.015 * H, f"{km:,} km", ha="center", va="bottom", fontsize=7, zorder=10)
        xa = x1 - 0.07 * W
        ax.annotate("", xy=(xa, y1 - 0.03 * H), xytext=(xa, y1 - 0.13 * H), arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2), zorder=10)
        ax.text(xa, y1 - 0.025 * H, "N", ha="center", va="bottom", fontsize=8, fontweight="bold", zorder=10)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)

    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.9), dpi=300)
    g.plot(column="rr_median", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=0.6, vcenter=1.0, vmax=1.6), ax=ax[0], linewidth=0, legend=True,
           legend_kwds={"shrink": 0.55, "label": "Posterior median RR", "orientation": "horizontal", "pad": 0.02})
    cuts = [0, 0.5, 0.8, 0.95, 1.0001]
    cols = ["#eff3ff", "#bdd7e7", "#6baed6", "#08306b"]
    g["cls"] = pd.cut(g.p_rr_gt1_10, cuts, labels=False, right=False)
    g.plot(column="cls", cmap=ListedColormap(cols), ax=ax[1], linewidth=0, categorical=False, vmin=0, vmax=3)
    for a, t in zip(ax, ["A", "B"]):
        states.boundary.plot(ax=a, linewidth=0.25, color="black"); a.set_axis_off(); a.set_title(t, loc="left", fontsize=9, fontweight="bold"); scale_and_north(a)
    ax[1].legend(handles=[Patch(facecolor=c, edgecolor="#555555", linewidth=0.4, label=l) for c, l in zip(cols, ["< 0.50", "0.50-0.79", "0.80-0.94", "≥ 0.95 (flagged)"])],
                 title="P(RR > 1.10)", frameon=False, loc="lower left", fontsize=6.5, title_fontsize=7)
    n_fl = int((g.p_rr_gt1_10 >= 0.95).sum())
    print("flagged on map:", n_fl, "| mapped:", len(g), "of", len(m))
    fig.tight_layout()
    fig.savefig(SUB / "Figure_1.png", dpi=300, bbox_inches="tight")
    fig.savefig(SUB / "Figure_1.pdf", bbox_inches="tight")
    plt.close(fig)
    print("figure 1 (map) ok")
else:
    print("MUNI_GEOJSON not set: figure 1 (map) skipped")
