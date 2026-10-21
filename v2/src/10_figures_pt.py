"""Portuguese-language manuscript figures (Figures 1 and 2) and supplementary Figure S1.

Same design as the English figures drawn in 05_results.py (figures), with Portuguese labels, decimal
commas, state names written in full (Figure S1), maps in the Brazilian polyconic projection (EPSG:5880),
and a scale bar (km) and north arrow on every map.

Reads:  supplement/municipal_estimates.csv (municipal results, distributed here) and the IBGE
        municipal and state meshes data/ref/br_municipios.geojson and br_states.geojson (not distributed).
Writes: outputs/figures_pt/figura1_rr_excedencia.{png,pdf}, figura2_obitos_excesso.{png,pdf},
        figura_s1_30_maiores_excessos.{png,pdf}

Run:
    python src/10_figures_pt.py
"""

from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
A1 = V2.parent
OUT = V2 / "outputs" / "figures_pt"
OUT.mkdir(parents=True, exist_ok=True)
CRS = 5880
UF = {"AC": "Acre", "AL": "Alagoas", "AP": "Amapá", "AM": "Amazonas", "BA": "Bahia", "CE": "Ceará",
      "DF": "Distrito Federal", "ES": "Espírito Santo", "GO": "Goiás", "MA": "Maranhão", "MT": "Mato Grosso",
      "MS": "Mato Grosso do Sul", "MG": "Minas Gerais", "PA": "Pará", "PB": "Paraíba", "PR": "Paraná",
      "PE": "Pernambuco", "PI": "Piauí", "RJ": "Rio de Janeiro", "RN": "Rio Grande do Norte",
      "RS": "Rio Grande do Sul", "RO": "Rondônia", "RR": "Roraima", "SC": "Santa Catarina",
      "SP": "São Paulo", "SE": "Sergipe", "TO": "Tocantins"}

m = pd.read_csv(V2 / "supplement" / "municipal_estimates.csv", dtype={"CODMUNRES": str})
munis = gpd.read_file(A1 / "data" / "ref" / "br_municipios.geojson").to_crs(CRS)
munis["CODMUNRES"] = munis["id"].astype(str).str[:6]
g = munis.merge(m, on="CODMUNRES", how="inner")
states = gpd.read_file(A1 / "data" / "ref" / "br_states.geojson").to_crs(CRS)
decimal_comma = FuncFormatter(lambda v, _: f"{v:.1f}".replace(".", ","))
thousands = lambda s: f"{s:,}".replace(",", ".")


def scale_and_north(ax, km=1000):
    """Scale bar (lower right) and north arrow (upper right), in projected metres."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); W, H = x1 - x0, y1 - y0
    L = km * 1000.0
    bx, by = x1 - 0.04 * W - L, y0 + 0.05 * H
    ax.plot([bx, bx + L], [by, by], color="black", lw=2.5, solid_capstyle="butt", zorder=10)
    for xx in (bx, bx + L):
        ax.plot([xx, xx], [by - 0.008 * H, by + 0.008 * H], color="black", lw=1.2, zorder=10)
    ax.text(bx + L / 2, by + 0.015 * H, f"{thousands(km)} km", ha="center", va="bottom", fontsize=8, zorder=10)
    xa = x1 - 0.07 * W
    ax.annotate("", xy=(xa, y1 - 0.03 * H), xytext=(xa, y1 - 0.13 * H),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5), zorder=10)
    ax.text(xa, y1 - 0.025 * H, "N", ha="center", va="bottom", fontsize=10, fontweight="bold", zorder=10)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)


fig, ax = plt.subplots(1, 2, figsize=(13, 6.5))
g.plot(column="rr_median", cmap="RdBu_r", vmin=0.6, vmax=1.6, ax=ax[0], linewidth=0, legend=True,
       legend_kwds={"shrink": 0.6, "label": "RR além da estrutura e da composição (mediana a posteriori)"})
states.boundary.plot(ax=ax[0], linewidth=0.3, color="black")
ax[0].set_axis_off(); ax[0].set_title("A. Risco relativo, mediana a posteriori")
g.plot(column="p_rr_gt1_10", cmap="viridis", vmin=0, vmax=1, ax=ax[1], linewidth=0, legend=True,
       legend_kwds={"shrink": 0.6, "label": "P(RR>1,10 | dados)"})
states.boundary.plot(ax=ax[1], linewidth=0.3, color="black")
ax[1].set_axis_off(); ax[1].set_title("B. Probabilidade de superar o esperado em mais de 10%")
for a in ax:
    scale_and_north(a)
for a in fig.axes:  # colour bars
    if a not in ax:
        a.yaxis.set_major_formatter(decimal_comma)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"figura1_rr_excedencia.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 9))
states.plot(ax=ax, color="#f2f2f2", edgecolor="black", linewidth=0.3)
ab = g[g["class"] == "above"].copy()
pts = ab.geometry.centroid
ax.scatter(pts.x, pts.y, s=np.clip(ab.excess_median, 1, None) / 6, color="#b2182b", alpha=0.6, edgecolor="black", linewidth=0.3)
for s_ in (50, 200, 800):
    ax.scatter([], [], s=s_ / 6, color="#b2182b", alpha=0.6, edgecolor="black", linewidth=0.3, label=f"{s_} óbitos")
ax.legend(loc="lower left", title="Óbitos em excesso (mediana a posteriori)", frameon=False)
ax.set_axis_off()
scale_and_north(ax)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"figura2_obitos_excesso.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)

top = m.sort_values("excess_median", ascending=False).head(30).iloc[::-1]
fig, ax = plt.subplots(figsize=(7, 8))
ax.errorbar(top.rr_median, range(len(top)), xerr=[top.rr_median - top.rr_lo95, top.rr_hi95 - top.rr_median],
            fmt="o", color="#b2182b", ecolor="gray", capsize=2)
ax.axvline(1.0, color="black", linewidth=0.8); ax.axvline(1.10, color="black", linewidth=0.8, linestyle="--")
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f"{a} ({UF[b]})  {thousands(int(c))} óbitos" for a, b, c in zip(top.municipio, top.state, top.deaths)], fontsize=7)
ax.xaxis.set_major_formatter(decimal_comma)
ax.set_xlabel("Risco relativo além da estrutura e da composição (mediana, ICr95%)")
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"figura_s1_30_maiores_excessos.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("municipalities above expectation:", len(ab), "| municipalities in the mesh:", len(g))
