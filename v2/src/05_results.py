"""Estimands, tables, figures and the numbers file for artigo 1 v2.

Everything the manuscript states comes from the posteriors written by 04_fit_models.py;
nothing is refitted here. The primary model is bym2_full.

Estimands (all with 95% equal-tailed credible intervals):
    RR_i           relative risk beyond structure and maternal composition, exp(u_i)
    exceedance     P(RR_i > 1), P(RR_i > 1.10), P(RR_i > 1.25), P(RR_i < 0.90)
    classification "above expectation" = P(RR_i > 1.10) >= 0.95; "below" = P(RR_i < 0.90) >= 0.95
    excess deaths  E_i (RR_i - 1) per municipality, where E_i = births_i exp(b0 + X_i beta);
                   summed as gross (positive part), net (signed), and above +10%, nationally,
                   by macro-region and by state
    heterogeneity  sigma_u (SD of log RR), rho (spatial share), sigma_hr (health region)
    detectability  minimum RR detectable with 80% power at a municipality's expected count
    negative control  sigma_u by cause group (malformations vs care-amenable groups)

Outputs:
    v2/outputs/tables/*.csv, v2/outputs/figures/*.png|pdf, v2/outputs/NUMBERS.json

Run:
    "../artigo 2/.venv/bin/python" v2/src/05_results.py
"""

import json
from pathlib import Path

import arviz as az
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from shapely.geometry import Point

matplotlib.use("Agg")
HERE = Path(__file__).resolve().parent
V2 = HERE.parent
A1 = V2.parent
DATA, MODEL, TAB, FIG = V2 / "data", V2 / "model", V2 / "outputs" / "tables", V2 / "outputs" / "figures"
TAB.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
PRIMARY = "bym2_rsr_full"
LABELS = {"idhm_renda": "IDHM income sub-index", "idhm_educ": "IDHM education sub-index", "logpop": "log population",
          "sh_age_lt20": "mothers < 20 y (share)", "sh_age_ge35": "mothers >= 35 y (share)",
          "sh_esc_le7": "mothers with <= 7 y schooling (share)", "sh_nonwhite": "non-white mothers (share)",
          "dist_nicu_100km": "distance to nearest NICU (per 100 km)", "has_nicu": "NICU beds in the municipality",
          "sh_born_outside": "births outside the municipality (share)", "sh_born_nicu_muni": "births in a NICU municipality (share)",
          "sh_cons_lt7": "< 7 antenatal visits (share)", "sh_preterm": "preterm births (share)", "sh_lbw": "low birth weight (share)",
          "sh_miss_esc": "births with missing schooling (share)", "sh_miss_race": "births with missing race/colour (share)",
          "sh_hospital": "hospital births (share)"}
REGIONS = ["North", "Northeast", "Central-West", "Southeast", "South"]
NUM = {}


def ci(a, axis=None):
    return np.quantile(a, [0.025, 0.975], axis=axis)


def post(name):
    return az.from_netcdf(MODEL / f"idata_{name}.nc")["posterior"]


def fit_info(name):
    return json.loads((MODEL / f"fit_{name}.json").read_text())


def hyper_table(df):
    """Hyperparameters and covariate effects (IRR per SD of the covariate) for every fitted model."""
    rows = []
    for p in sorted(MODEL.glob("idata_*.nc")):
        name = p.stem.replace("idata_", "")
        if name.startswith("period"):
            continue  # period-specific fits are summarised by 07_period_persistence.py
        ps = post(name)
        info = fit_info(name)
        X = info["spec"]["X"]
        sd = {c: float(df[c].std()) for c in X}
        beta = ps["beta"].stack(s=("chain", "draw")).values
        for k, c in enumerate(X):
            irr = np.exp(beta[k])
            rows.append({"model": name, "term": LABELS.get(c, c), "covariate": c, "unit": f"per SD ({sd[c]:.3g})",
                         "irr": float(np.median(irr)), "lo95": float(ci(irr)[0]), "hi95": float(ci(irr)[1]),
                         "p_gt1": float((irr > 1).mean())})
        for h in ["sigma_u", "rho", "sigma_hr"]:
            if h in ps:
                v = ps[h].values.ravel()
                rows.append({"model": name, "term": h, "covariate": h, "unit": "",
                             "irr": float(np.median(v)), "lo95": float(ci(v)[0]), "hi95": float(ci(v)[1]), "p_gt1": np.nan})
        rows.append({"model": name, "term": "diagnostics", "covariate": "", "unit": "",
                     "irr": info["divergences"], "lo95": info["max_rhat_u"], "hi95": info["min_ess_bulk_u"], "p_gt1": np.nan})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "model_hyperparameters_all.csv", index=False)
    return t


def primary_estimands(df):
    ps = post(PRIMARY)
    u = ps["u"].stack(s=("chain", "draw")).values            # n x S
    b0 = ps["b0"].values.ravel()                               # S
    beta = ps["beta"].stack(s=("chain", "draw")).values        # k x S
    X = df[fit_info(PRIMARY)["spec"]["X"]].values.astype(float)
    Xz = (X - X.mean(0)) / X.std(0)
    births = df.live_births.values.astype(float)
    y = df.deaths_total.values.astype(float)
    E = births[:, None] * np.exp(b0[None, :] + Xz @ beta)      # n x S expected from covariates
    rr = np.exp(u)
    excess = E * (rr - 1.0)                                    # signed, n x S
    sigma = ps["sigma_u"].values.ravel()
    rho = ps["rho"].values.ravel()

    m = pd.DataFrame({
        "CODMUNRES": df.CODMUNRES, "municipio": df.municipio, "state": df.state, "REG": df.REG,
        "deaths": y.astype(int), "births": births.astype(int), "nmr": y / births * 1000,
        "expected_median": np.median(E, 1), "expected_lo95": ci(E, 1)[0], "expected_hi95": ci(E, 1)[1],
        "rr_median": np.median(rr, 1), "rr_lo95": ci(rr, 1)[0], "rr_hi95": ci(rr, 1)[1],
        "p_rr_gt1": (rr > 1).mean(1), "p_rr_gt1_10": (rr > 1.10).mean(1), "p_rr_gt1_25": (rr > 1.25).mean(1),
        "p_rr_lt0_90": (rr < 0.90).mean(1),
        "excess_median": np.median(excess, 1), "excess_lo95": ci(excess, 1)[0], "excess_hi95": ci(excess, 1)[1],
        "v1_smr": df.v1_smr, "v1_sig_excess": df.v1_sig_excess, "sh_hospital": df.sh_hospital,
        "sh_miss_esc": df.sh_miss_esc, "sh_miss_race": df.sh_miss_race,
    })
    m["class"] = np.select([m.p_rr_gt1_10 >= 0.95, m.p_rr_lt0_90 >= 0.95], ["above", "below"], "undetermined")
    m.sort_values("rr_median", ascending=False).to_csv(TAB / "municipal_estimates.csv", index=False)

    # National and regional sums with full posterior propagation.
    def sums(mask, label):
        g = excess[mask].clip(min=0).sum(0); n_ = excess[mask].sum(0)
        a10 = (E[mask] * (rr[mask] - 1.10)).clip(min=0).sum(0)
        above = (m["class"].values == "above") & mask
        e_above = excess[above].sum(0) if above.any() else np.zeros(excess.shape[1])
        return {"unit": label, "municipalities": int(mask.sum()), "deaths": int(y[mask].sum()),
                "expected_median": float(np.median(E[mask].sum(0))),
                "gross_median": float(np.median(g)), "gross_lo95": float(ci(g)[0]), "gross_hi95": float(ci(g)[1]),
                "net_median": float(np.median(n_)), "net_lo95": float(ci(n_)[0]), "net_hi95": float(ci(n_)[1]),
                "above10_median": float(np.median(a10)), "above10_lo95": float(ci(a10)[0]), "above10_hi95": float(ci(a10)[1]),
                "n_above": int(above.sum()), "n_below": int(((m["class"].values == "below") & mask).sum()),
                "excess_in_above_median": float(np.median(e_above)), "excess_in_above_lo95": float(ci(e_above)[0]),
                "excess_in_above_hi95": float(ci(e_above)[1]),
                "gross_share_of_deaths": float(np.median(g) / y[mask].sum() * 100)}
    rows = [sums(np.ones(len(df), bool), "Brazil")]
    rows += [sums(df.REG.values == r, r) for r in REGIONS]
    rows += [sums(df.state.values == s, s) for s in sorted(df.state.unique())]
    agg = pd.DataFrame(rows)
    agg.to_csv(TAB / "excess_by_unit.csv", index=False)

    # What a world with no residual heterogeneity would give the same summaries: the gross
    # positive part under pure Poisson noise around E, to show the floor the v1 ladder sat on.
    rng = np.random.default_rng(1)
    S = min(500, E.shape[1])
    idx = rng.choice(E.shape[1], S, replace=False)
    ysim = rng.poisson(E[:, idx])
    floor = (ysim - E[:, idx]).clip(min=0).sum(0)
    NUM["poisson_floor_gross_all"] = {"median": float(np.median(floor)), "lo95": float(ci(floor)[0]), "hi95": float(ci(floor)[1])}
    rel = y >= 50
    floor_rel = (ysim[rel] - E[rel][:, idx]).clip(min=0).sum(0)
    NUM["poisson_floor_gross_ge50"] = {"median": float(np.median(floor_rel)), "lo95": float(ci(floor_rel)[0]), "hi95": float(ci(floor_rel)[1])}

    # Top municipalities.
    top = m.sort_values("rr_median", ascending=False).head(20)
    top.to_csv(TAB / "top20_above.csv", index=False)
    topx = m[m["class"] == "above"].sort_values("excess_median", ascending=False).head(20)
    topx.to_csv(TAB / "top20_excess_deaths.csv", index=False)

    sd_real = u.std(0)
    NUM.update({
        "n_municipalities": int(len(df)), "deaths_modelled": int(y.sum()), "births_modelled": int(births.sum()),
        "sigma_u": {"median": float(np.median(sigma)), "lo95": float(ci(sigma)[0]), "hi95": float(ci(sigma)[1])},
        "rho": {"median": float(np.median(rho)), "lo95": float(ci(rho)[0]), "hi95": float(ci(rho)[1])},
        "u_sd_realised": {"median": float(np.median(sd_real)), "lo95": float(ci(sd_real)[0]), "hi95": float(ci(sd_real)[1])},
        "rr_sd_implied_pct": float((np.exp(np.median(sd_real)) - 1) * 100),
        "n_above_10": int((m["class"] == "above").sum()), "n_below_10": int((m["class"] == "below").sum()),
        "n_p_gt1_25_ge95": int((m.p_rr_gt1_25 >= 0.95).sum()), "n_p_gt1_ge95": int((m.p_rr_gt1 >= 0.95).sum()),
        "n_p_gt1_10_ge80": int((m.p_rr_gt1_10 >= 0.80).sum()),
        "above_by_region": m[m["class"] == "above"].REG.value_counts().to_dict(),
        "below_by_region": m[m["class"] == "below"].REG.value_counts().to_dict(),
        "v1_flagged": int(df.v1_sig_excess.sum()),
        "v1_flagged_now_above": int(((m["class"] == "above") & df.v1_sig_excess).sum()),
        "v1_flagged_now_p10_lt50": int(((m.p_rr_gt1_10 < 0.5) & df.v1_sig_excess).sum()),
        "above_not_flagged_v1": int(((m["class"] == "above") & ~df.v1_sig_excess.astype(bool)).sum()),
        "v1_reliable_flagged": int((df.v1_sig_excess & rel).sum()),
        "v1_reliable_flagged_now_above": int(((m["class"] == "above") & df.v1_sig_excess & rel).sum()),
        "v1_reliable_flagged_now_p10_lt50": int(((m.p_rr_gt1_10 < 0.5) & df.v1_sig_excess & rel).sum()),
        "v1_reliable_flagged_now_p10_lt80": int(((m.p_rr_gt1_10 < 0.8) & df.v1_sig_excess & rel).sum()),
        "above_small": int(((m["class"] == "above") & ~rel).sum()),
        "above_deaths_total": int(y[(m["class"] == "above").values].sum()),
        "above_excess_share_of_their_deaths": float(rows[0]["excess_in_above_median"] / y[(m["class"] == "above").values].sum() * 100),
        "below_deaths_total": int(y[(m["class"] == "below").values].sum()),
        "above_by_state": m[m["class"] == "above"].state.value_counts().to_dict(),
        "below_by_state": m[m["class"] == "below"].state.value_counts().to_dict(),
        "registration_above_vs_other_median": {c: [float(m.loc[m["class"] == "above", c].median()), float(m.loc[m["class"] != "above", c].median())]
                                               for c in ["sh_hospital", "sh_miss_esc", "sh_miss_race"]},
        "above_low_hospital": m.loc[(m["class"] == "above") & (m.sh_hospital < 0.90), ["municipio", "state", "deaths", "rr_median", "sh_hospital", "sh_miss_esc", "sh_miss_race"]]
                               .sort_values("rr_median", ascending=False).to_dict("records"),
        "national_hospital_share": float((df.sh_hospital * df.births).sum() / df.births.sum()),
        "rr_median_iqr_all": [float(m.rr_median.quantile(q)) for q in (.25, .5, .75)],
        "rr_median_range_above": [float(m[m["class"] == "above"].rr_median.min()), float(m[m["class"] == "above"].rr_median.max())],
        "spearman_rr_v1smr_ge50": float(stats.spearmanr(m.rr_median[rel], m.v1_smr[rel])[0]),
        "spearman_rr_hospital": float(stats.spearmanr(m.rr_median, m.sh_hospital)[0]),
        "spearman_rr_miss_race": float(stats.spearmanr(m.rr_median, m.sh_miss_race)[0]),
        "brazil": rows[0], "regions": {r["unit"]: r for r in rows[1:6]},
    })
    return m, E, rr, agg


def expected_and_rr(df, name):
    """Expected deaths from the covariates and the residual relative risk, draw by draw."""
    ps = post(name)
    u = ps["u"].stack(s=("chain", "draw")).values.astype(float)
    b0 = ps["b0"].values.ravel()
    beta = ps["beta"].stack(s=("chain", "draw")).values
    X = df[fit_info(name)["spec"]["X"]].values.astype(float)
    Xz = (X - X.mean(0)) / X.std(0)
    E = df.live_births.values.astype(float)[:, None] * np.exp(b0[None, :] + Xz @ beta)
    return u, E, np.exp(u)


def waic(y, mu):
    """Widely applicable information criterion from pointwise Poisson log-likelihoods."""
    ll = stats.poisson.logpmf(y[:, None], mu)
    mx = ll.max(1, keepdims=True)
    lppd = (np.log(np.exp(ll - mx).mean(1)) + mx.ravel())
    p = ll.var(1)
    return float((lppd - p).sum()), float(p.sum())


def model_comparison(df, m):
    """How the heterogeneity, the fit and the map move across covariate sets and priors.

    The heterogeneity reported is the realised spread of u across municipalities within a
    draw (its SD, averaged over draws), not the BYM2 scale parameter sigma_u: under a scaled
    ICAR the two differ, and the realised spread is what the map shows.
    """
    rows = []
    y = df.deaths_total.values.astype(float)
    ref = post(PRIMARY)["u"].stack(s=("chain", "draw")).values.mean(1)
    for name in ["iid_struct", "iid_full", "bym2_rsr_struct", PRIMARY, "bym2_rsr_full_hr", "bym2_struct", "bym2_full",
                 "bym2_full_hr", "iid_full_access", "iid_full_intermediate", "iid_full_completeness"]:
        p = MODEL / f"idata_{name}.nc"
        if not p.exists():
            continue
        ps = post(name)
        s = ps["sigma_u"].values.ravel()
        um, E, rr = expected_and_rr(df, name)
        sd_real = um.std(0)
        elpd, p_eff = waic(y, E * rr)
        gross = (E * (rr - 1)).clip(min=0).sum(0)
        row = {"model": name, "sigma_u_median": float(np.median(s)), "sigma_u_lo95": float(ci(s)[0]), "sigma_u_hi95": float(ci(s)[1]),
               "u_sd_realised_median": float(np.median(sd_real)), "u_sd_realised_lo95": float(ci(sd_real)[0]),
               "u_sd_realised_hi95": float(ci(sd_real)[1]), "elpd_waic": elpd, "p_waic": p_eff,
               "gross_median": float(np.median(gross)), "gross_lo95": float(ci(gross)[0]), "gross_hi95": float(ci(gross)[1])}
        if "rho" in ps:
            r = ps["rho"].values.ravel(); row.update(rho_median=float(np.median(r)), rho_lo95=float(ci(r)[0]), rho_hi95=float(ci(r)[1]))
        if "sigma_hr" in ps:
            # Share of the realised residual variance that sits at the health-region level.
            h = ps["sigma_hr"].values.ravel()
            icc = h ** 2 / sd_real ** 2
            row.update(sigma_hr_median=float(np.median(h)), sigma_hr_lo95=float(ci(h)[0]), sigma_hr_hi95=float(ci(h)[1]),
                       icc_hr_median=float(np.median(icc)), icc_hr_lo95=float(ci(icc)[0]), icc_hr_hi95=float(ci(icc)[1]))
        row.update(n_above_10=int(((rr > 1.10).mean(1) >= 0.95).sum()), n_below_10=int(((rr < 0.90).mean(1) >= 0.95).sum()),
                   spearman_u_vs_primary=float(stats.spearmanr(um.mean(1), ref)[0]))
        rows.append(row)
    t = pd.DataFrame(rows)
    t["delta_elpd_vs_primary"] = t.elpd_waic - float(t.loc[t.model == PRIMARY, "elpd_waic"].iloc[0])
    t.to_csv(TAB / "model_comparison.csv", index=False)
    NUM["model_comparison"] = t.set_index("model").to_dict("index")
    return t


def cause_control(df, m):
    rows = []
    for c in ["prenatal", "delivery", "newborn", "malformation"]:
        p = MODEL / f"idata_cause_{c}.nc"
        if not p.exists():
            continue
        ps = post(f"cause_{c}")
        s = ps["sigma_u"].values.ravel()
        rr = np.exp(ps["u"].stack(s=("chain", "draw")).values)
        above = (m["class"] == "above").values
        rows.append({"cause": c, "deaths": int(df[c].sum()), "sigma_u_median": float(np.median(s)),
                     "sigma_u_lo95": float(ci(s)[0]), "sigma_u_hi95": float(ci(s)[1]),
                     "rr_median_in_above_set": float(np.median(np.median(rr[above], 1))),
                     "rr_median_in_others": float(np.median(np.median(rr[~above], 1))),
                     "n_above_10": int(((rr > 1.10).mean(1) >= 0.95).sum())})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "cause_negative_control.csv", index=False)
    NUM["cause_control"] = t.set_index("cause").to_dict("index")
    return t


def detectability(df, E):
    sigma = float(NUM["u_sd_realised"]["median"])
    z = stats.norm.ppf(0.95) + stats.norm.ppf(0.80)
    Em = np.median(E, 1)
    rows = []
    for e in [10, 20, 50, 100, 200, 500, 1000, 2000, 4000]:
        mde_p = 1 + z * np.sqrt(e) / e
        mde_h = 1 + z * np.sqrt(e + (np.exp(sigma ** 2) - 1) * e ** 2) / e
        rows.append({"expected_deaths": e, "mde_rr_poisson": mde_p, "mde_rr_with_heterogeneity": mde_h,
                     "municipalities_at_or_above": int((Em >= e).sum()), "deaths_share_at_or_above": float(df.deaths_total[Em >= e].sum() / df.deaths_total.sum())})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "detectability.csv", index=False)
    NUM["detectability"] = t.to_dict("records")
    NUM["mde_asymptote"] = float(1 + z * np.sqrt(np.exp(sigma ** 2) - 1))
    return t


def ppc(df, E, rr):
    rng = np.random.default_rng(2)
    S = min(1000, E.shape[1])
    idx = rng.choice(E.shape[1], S, replace=False)
    mu = E[:, idx] * rr[:, idx]
    ysim = rng.poisson(mu)
    y = df.deaths_total.values
    lo, hi = np.quantile(ysim, [0.025, 0.975], axis=1)
    out = {"share_outside_95pp": float(((y < lo) | (y > hi)).mean()),
           "share_zero_obs": float((y == 0).mean()), "share_zero_pp": float((ysim == 0).mean()),
           "max_obs": int(y.max()), "max_pp_median": float(np.median(ysim.max(0))),
           "var_mean_ratio_obs": float(y.var() / y.mean()), "var_mean_ratio_pp_median": float(np.median(ysim.var(0) / ysim.mean(0)))}
    NUM["ppc"] = out
    return out


def figures(df, m):
    munis = gpd.read_file(A1 / "data" / "ref" / "br_municipios.geojson")
    munis["CODMUNRES"] = munis["id"].astype(str).str[:6]
    g = munis.merge(m, on="CODMUNRES", how="inner")
    states = gpd.read_file(A1 / "data" / "ref" / "br_states.geojson")

    fig, ax = plt.subplots(1, 2, figsize=(13, 6.5))
    g.plot(column="rr_median", cmap="RdBu_r", vmin=0.6, vmax=1.6, ax=ax[0], linewidth=0, legend=True,
           legend_kwds={"shrink": 0.6, "label": "Relative risk beyond structure and composition (posterior median)"})
    states.boundary.plot(ax=ax[0], linewidth=0.3, color="black")
    ax[0].set_axis_off(); ax[0].set_title("A. Relative risk, posterior median")
    g.plot(column="p_rr_gt1_10", cmap="viridis", vmin=0, vmax=1, ax=ax[1], linewidth=0, legend=True,
           legend_kwds={"shrink": 0.6, "label": "P(RR > 1.10 | data)"})
    states.boundary.plot(ax=ax[1], linewidth=0.3, color="black")
    ax[1].set_axis_off(); ax[1].set_title("B. Probability of exceeding the expectation by more than 10%")
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figure1_rr_exceedance.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 9))
    states.plot(ax=ax, color="#f2f2f2", edgecolor="black", linewidth=0.3)
    ab = g[g["class"] == "above"].copy()
    pts = ab.to_crs(5880).geometry.centroid.to_crs(4326)
    ax.scatter(pts.x, pts.y, s=np.clip(ab.excess_median, 1, None) / 6, color="#b2182b", alpha=0.6, edgecolor="black", linewidth=0.3)
    for s_ in (50, 200, 800):
        ax.scatter([], [], s=s_ / 6, color="#b2182b", alpha=0.6, edgecolor="black", linewidth=0.3, label=f"{s_} deaths")
    ax.legend(loc="lower left", title="Excess deaths (posterior median)", frameon=False)
    ax.set_axis_off(); ax.set_title("Municipalities above expectation with >= 95% probability (RR > 1.10)")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figure2_excess_deaths.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    top = m.sort_values("excess_median", ascending=False).head(30).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 8))
    ax.errorbar(top.rr_median, range(len(top)), xerr=[top.rr_median - top.rr_lo95, top.rr_hi95 - top.rr_median],
                fmt="o", color="#b2182b", ecolor="gray", capsize=2)
    ax.axvline(1.0, color="black", linewidth=0.8); ax.axvline(1.10, color="black", linewidth=0.8, linestyle="--")
    ax.set_yticks(range(len(top))); ax.set_yticklabels([f"{a} ({b})  {int(c):,} deaths" for a, b, c in zip(top.municipio, top.state, top.deaths)], fontsize=7)
    ax.set_xlabel("Relative risk beyond structure and composition (median, 95% CrI)")
    ax.set_title("The 30 municipalities with the largest excess deaths")
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figure3_caterpillar.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main():
    df = pd.read_csv(DATA / "analysis_muni.csv", dtype={"CODMUNRES": str, "regsaud_id": str})
    hyper = hyper_table(df)
    m, E, rr, agg = primary_estimands(df)
    model_comparison(df, m)
    cause_control(df, m)
    detectability(df, E)
    ppc(df, E, rr)
    figures(df, m)
    NUM["primary_hyper"] = hyper[hyper.model == PRIMARY].drop(columns="model").to_dict("records")
    NUM["fit_primary"] = {k: v for k, v in fit_info(PRIMARY).items() if k != "hyper"}
    (V2 / "outputs" / "NUMBERS.json").write_text(json.dumps(NUM, indent=2, default=float))
    print(json.dumps({k: NUM[k] for k in ["sigma_u", "rho", "n_above_10", "n_below_10", "n_p_gt1_25_ge95", "brazil",
                                          "v1_flagged", "v1_flagged_now_above", "above_not_flagged_v1", "ppc", "mde_asymptote"]},
                     indent=1, default=float))
    print(agg[agg.unit.isin(["Brazil"] + REGIONS)][["unit", "deaths", "gross_median", "gross_lo95", "gross_hi95", "net_median",
                                                       "net_lo95", "net_hi95", "n_above", "n_below", "excess_in_above_median"]].round(0).to_string(index=False))


if __name__ == "__main__":
    main()
