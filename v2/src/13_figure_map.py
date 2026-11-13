"""Figure 3 (English): residual relative risk map, derived from v2/src/10_figures_pt.py (same design, English labels,
decimal points, scale bar and north arrow; Calibri 12 pt). Reads v2/outputs/tables/municipal_estimates.csv and the IBGE
meshes in data/ref (not distributed)."""
import glob
from pathlib import Path
import geopandas as gpd, matplotlib, pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib import font_manager as fm
for _f in glob.glob('/Applications/Microsoft Excel.app/Contents/Resources/DFonts/Calibri*.ttf'):
    fm.fontManager.addfont(_f)
names = {f.name for f in fm.fontManager.ttflist}
plt.rcParams['font.family'] = next((n for n in ('Calibri','Carlito','Arial','DejaVu Sans') if n in names), 'DejaVu Sans')
plt.rcParams['font.size'] = 12
print('fonte:', plt.rcParams['font.family'])
HERE = Path(__file__).resolve().parent; V2 = HERE.parent; A1 = V2.parent; OUT = V2 / "outputs" / "figures"; OUT.mkdir(parents=True, exist_ok=True)
CRS = 5880
m = pd.read_csv(V2 / "outputs" / "tables" / "municipal_estimates.csv", dtype={"CODMUNRES": str})
munis = gpd.read_file(A1 / "data" / "ref" / "br_municipios.geojson").to_crs(CRS)
munis["CODMUNRES"] = munis["id"].astype(str).str[:6]
g = munis.merge(m, on="CODMUNRES", how="inner")
states = gpd.read_file(A1 / "data" / "ref" / "br_states.geojson").to_crs(CRS)
print("RR above 1.6 (shown in the darkest red):", int((g.rr_median > 1.6).sum()), "| max RR:", round(float(g.rr_median.max()), 3))
print("municipalities in the estimates:", len(m), "| in the mesh and estimates:", len(g))

def scale_and_north(ax, km=1000):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); W, H = x1 - x0, y1 - y0
    L = km * 1000.0
    bx, by = x1 - 0.04 * W - L, y0 + 0.05 * H
    ax.plot([bx, bx + L], [by, by], color="black", lw=2.5, solid_capstyle="butt", zorder=10)
    for xx in (bx, bx + L):
        ax.plot([xx, xx], [by - 0.008 * H, by + 0.008 * H], color="black", lw=1.2, zorder=10)
    ax.text(bx + L / 2, by + 0.015 * H, f"{km:,} km", ha="center", va="bottom", fontsize=12, zorder=10)
    xa = x1 - 0.07 * W
    ax.annotate("", xy=(xa, y1 - 0.03 * H), xytext=(xa, y1 - 0.13 * H),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5), zorder=10)
    ax.text(xa, y1 - 0.025 * H, "N", ha="center", va="bottom", fontsize=12, fontweight="bold", zorder=10)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)

fig, ax = plt.subplots(1, 2, figsize=(13, 6.5))
g.plot(column="rr_median", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=0.6, vcenter=1.0, vmax=1.6), ax=ax[0], linewidth=0, legend=True,
       legend_kwds={"shrink": 0.6, "label": "Relative risk beyond structure and composition (posterior median)"})
states.boundary.plot(ax=ax[0], linewidth=0.3, color="black"); ax[0].set_axis_off(); ax[0].set_title("A", loc="left", fontsize=14, fontweight="bold")
g.plot(column="p_rr_gt1_10", cmap="viridis", vmin=0, vmax=1, ax=ax[1], linewidth=0, legend=True,
       legend_kwds={"shrink": 0.6, "label": "P(RR > 1.10 | data)"})
states.boundary.plot(ax=ax[1], linewidth=0.3, color="black"); ax[1].set_axis_off(); ax[1].set_title("B", loc="left", fontsize=14, fontweight="bold")
for a in ax: scale_and_north(a)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"figure3_rr_map.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
