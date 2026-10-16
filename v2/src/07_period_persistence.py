"""Persistence of the residual relative risk between 2014-2019 and 2020-2024.

The primary model (BYM2, restricted spatial regression, structure + maternal composition) is fitted
separately to each period, with the same 5,565 municipalities, the same graph, the same priors and
the same sampler settings. Period-specific inputs: deaths (year of death), live births (year of birth),
resident population (mean over the period) and maternal composition (births of the period). The HDI
(2010) is the same in both periods.

Questions: do the same municipalities lie above expectation in both periods, and is the residual
ordering of municipalities stable between periods?

Inputs:  data/analysis_muni.csv, data/adjacency_pairs.csv, data/bym2_scaling.json,
         data/maternal_composition_muni_year.csv, data/processed/nascidos_muni_ano.csv,
         data/processed/populacao_muni_ano.csv, data/deaths_muni_year.csv
         (built once from the individual-level SIM neonatal extract in data/raw/ if absent)
Outputs: outputs/tables/period_estimates.csv, outputs/period_persistence.json, model/idata_period*.nc

Run:
    python src/07_period_persistence.py
"""

import importlib.util
import json
import os as _os
import time
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import xarray as xr
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = Path(_os.environ.get("NNM1_ROOT", HERE.parent))  # repository root; override with NNM1_ROOT


def first_existing(*candidates):
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


PROC = first_existing(ROOT / "data" / "processed", ROOT.parent / "data" / "processed")
RAW = first_existing(ROOT / "data" / "raw", ROOT.parent / "data" / "raw")
DATA, MODEL = ROOT / "data", ROOT / "model"
TAB = ROOT / "outputs" / "tables"
MODEL.mkdir(exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

_spec = importlib.util.spec_from_file_location("fit_models", HERE / "04_fit_models.py")
fm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fm)

PERIODS = {1: (2014, 2019), 2: (2020, 2024)}
X = fm.STRUCT + fm.COMPOSITION


def deaths_by_year():
    path = DATA / "deaths_muni_year.csv"
    if path.exists():
        return pd.read_csv(path, dtype={"CODMUNRES": str})
    raw = pd.read_csv(RAW / "coorte_neonatal_2014_2024.csv", sep=";", dtype=str, usecols=["CODMUNRES", "DTOBITO"])
    raw["CODMUNRES"] = raw.CODMUNRES.str[:6]
    raw = raw[~raw.CODMUNRES.str.endswith("0000")]
    raw["year"] = raw.DTOBITO.str[-4:].astype(int)
    g = raw.groupby(["CODMUNRES", "year"]).size().rename("deaths").reset_index()
    g.to_csv(path, index=False)
    return g


def period_frame(base, period, dy, births, pop, comp):
    y0, y1 = PERIODS[period]
    df = base[["CODMUNRES", "municipio", "state", "REG", "regsaud_id", "idhm_renda", "idhm_educ"]].copy()
    d = dy[(dy.year >= y0) & (dy.year <= y1)].groupby("CODMUNRES").deaths.sum()
    b = births[(births.ANO >= y0) & (births.ANO <= y1)].groupby("CODMUNRES").nascimentos.sum()
    p = pop[(pop.ANO >= y0) & (pop.ANO <= y1)].groupby("CODMUNRES").populacao.mean()
    c = comp[(comp.year >= y0) & (comp.year <= y1)].drop(columns="year").groupby("cod6").sum()
    df["deaths_total"] = df.CODMUNRES.map(d).fillna(0).astype(int)
    df["live_births"] = df.CODMUNRES.map(b).fillna(0).astype(int)
    df["logpop"] = np.log(df.CODMUNRES.map(p))
    for new, num, den in [("sh_age_lt20", "age_lt20", "age_known"), ("sh_age_ge35", "age_ge35", "age_known"),
                          ("sh_esc_le7", "esc_le7", "esc_known"), ("sh_nonwhite", "race_nonwhite", "race_known")]:
        n_, d_ = df.CODMUNRES.map(c[num]), df.CODMUNRES.map(c[den])
        df[new] = np.where(d_ > 0, n_ / d_.replace(0, np.nan), np.nan)
        df[new] = df[new].fillna(df[new].median())
    assert (df.live_births > 0).all(), "a municipality has no births in the period"
    return df


def fit_period(df, W, scaling, period, draws=1000, tune=2000, chains=4):
    t0 = time.time()
    with fm.build(df, W, scaling, fm.SPECS["bym2_rsr_full"]):
        raw = pm.sample(draws=draws, tune=tune, chains=chains, cores=chains, random_seed=fm.SEED + period,
                        target_accept=0.9, progressbar=False)
    ds = raw["posterior"].to_dataset()[["b0", "beta", "sigma_u", "rho", "u"]]
    ds["u"] = ds["u"].astype("float32")
    idata = xr.DataTree.from_dict({"posterior": ds, "sample_stats": raw["sample_stats"].to_dataset()})
    idata.to_netcdf(MODEL / f"idata_period{period}.nc", engine="h5netcdf")
    info = {"divergences": int(idata["sample_stats"]["diverging"].sum()),
            "max_rhat_u": float(az.rhat(idata, var_names=["u"])["u"].max()),
            "elapsed_min": round((time.time() - t0) / 60, 1)}
    print(f"period {period}: {info}", flush=True)
    return idata, info


def draws_of(idata, df):
    post = idata["posterior"]
    u = post["u"].stack(s=("chain", "draw")).values
    b0 = post["b0"].stack(s=("chain", "draw")).values
    beta = post["beta"].stack(s=("chain", "draw")).values
    Xr = df[X].values.astype(float)
    Xz = (Xr - Xr.mean(0)) / Xr.std(0)
    E = df.live_births.values[:, None] * np.exp(b0[None, :] + Xz @ beta)
    return u, E, post["sigma_u"].values.ravel(), post["rho"].values.ravel()


def ci(a, axis=None):
    return np.quantile(a, [0.025, 0.975], axis=axis)


def main():
    base = pd.read_csv(DATA / "analysis_muni.csv", dtype={"CODMUNRES": str, "regsaud_id": str})
    pairs = pd.read_csv(DATA / "adjacency_pairs.csv")
    scaling = json.loads((DATA / "bym2_scaling.json").read_text())["scaling_factor"]
    n = len(base)
    W = np.zeros((n, n))
    W[pairs.i.values, pairs.j.values] = 1.0
    W[pairs.j.values, pairs.i.values] = 1.0
    dy = deaths_by_year()
    births = pd.read_csv(PROC / "nascidos_muni_ano.csv", dtype={"CODMUNRES": str})
    births["CODMUNRES"] = births.CODMUNRES.str[:6]
    pop = pd.read_csv(PROC / "populacao_muni_ano.csv", dtype={"CODMUNRES": str})
    comp = pd.read_csv(DATA / "maternal_composition_muni_year.csv", dtype={"cod6": str})

    res, frames = {}, {}
    for p in PERIODS:
        frames[p] = period_frame(base, p, dy, births, pop, comp)
        path = MODEL / f"idata_period{p}.nc"
        if path.exists():
            idata = az.from_netcdf(path)
            prev = ROOT / "outputs" / "period_persistence.json"
            info = json.loads(prev.read_text()).get("fit", {}).get(str(p), {"cached": True}) if prev.exists() else {"cached": True}
        else:
            idata, info = fit_period(frames[p], W, scaling, p)
        u, E, sig, rho = draws_of(idata, frames[p])
        res[p] = dict(u=u, E=E, sig=sig, rho=rho, info=info)

    out = pd.DataFrame({"CODMUNRES": base.CODMUNRES, "municipio": base.municipio, "state": base.state, "REG": base.REG})
    flag, rr_med, p110 = {}, {}, {}
    for p in PERIODS:
        rr = np.exp(res[p]["u"])
        out[f"deaths_p{p}"] = frames[p].deaths_total
        out[f"births_p{p}"] = frames[p].live_births
        out[f"expected_p{p}"] = np.median(res[p]["E"], 1)
        out[f"rr_median_p{p}"] = np.median(rr, 1)
        out[f"rr_lo95_p{p}"], out[f"rr_hi95_p{p}"] = ci(rr, 1)
        out[f"p_rr_gt1_10_p{p}"] = (rr > 1.10).mean(1)
        out[f"p_rr_lt0_90_p{p}"] = (rr < 0.90).mean(1)
        out[f"class_p{p}"] = np.select([out[f"p_rr_gt1_10_p{p}"] >= 0.95, out[f"p_rr_lt0_90_p{p}"] >= 0.95], ["above", "below"], "undetermined")
        flag[p] = (out[f"class_p{p}"] == "above").values
        rr_med[p] = out[f"rr_median_p{p}"].values
        p110[p] = out[f"p_rr_gt1_10_p{p}"].values
    out.to_csv(TAB / "period_estimates.csv", index=False)

    pooled = pd.read_csv(TAB / "municipal_estimates.csv", dtype={"CODMUNRES": str}).set_index("CODMUNRES")
    pooled_rr = pooled.rr_median.reindex(out.CODMUNRES).values
    pooled_flag = (pooled["class"].reindex(out.CODMUNRES) == "above").values
    both, only1, only2 = flag[1] & flag[2], flag[1] & ~flag[2], ~flag[1] & flag[2]
    n1, n2, N = int(flag[1].sum()), int(flag[2].sum()), n
    exp_both = n1 * n2 / N
    hyp_p = float(stats.hypergeom.sf(int(both.sum()) - 1, N, n1, n2))
    # size-stratified permutation benchmark: flags of period 2 are shuffled within strata of pooled expected deaths
    exp_pooled = pooled.expected_median.reindex(out.CODMUNRES).values
    strata = np.digitize(exp_pooled, [20, 50, 100, 200, 500])
    rng = np.random.default_rng(20261016)
    perm = np.empty(10000, dtype=int)
    for b in range(perm.size):
        f2 = np.zeros(N, dtype=bool)
        for sidx in np.unique(strata):
            idx = np.where(strata == sidx)[0]
            k2 = int(flag[2][idx].sum())
            if k2:
                f2[rng.choice(idx, k2, replace=False)] = True
        perm[b] = int((flag[1] & f2).sum())
    big = (out.deaths_p1 >= 25) & (out.deaths_p2 >= 25)
    big20 = (out.expected_p1 >= 20) & (out.expected_p2 >= 20)
    sp = lambda a, b, m=slice(None): float(stats.spearmanr(a[m], b[m])[0])
    S = {
        "periods": {"1": "2014-2019", "2": "2020-2024"},
        "deaths": {p: int(frames[p].deaths_total.sum()) for p in PERIODS},
        "births": {p: int(frames[p].live_births.sum()) for p in PERIODS},
        "fit": {p: res[p]["info"] for p in PERIODS},
        "sd_log_rr": {p: {"median": float(np.median(res[p]["u"].std(0))),
                          "lo95": float(ci(res[p]["u"].std(0))[0]), "hi95": float(ci(res[p]["u"].std(0))[1])} for p in PERIODS},
        "rho": {p: float(np.median(res[p]["rho"])) for p in PERIODS},
        "n_above": {p: int(flag[p].sum()) for p in PERIODS},
        "n_below": {p: int((out[f"class_p{p}"] == "below").sum()) for p in PERIODS},
        "above_both": int(both.sum()), "above_p1_only": int(only1.sum()), "above_p2_only": int(only2.sum()),
        "above_both_expected_if_independent": float(exp_both), "above_both_hypergeom_p": hyp_p,
        "above_both_stratified_permutation_mean": float(perm.mean()), "above_both_stratified_permutation_max": int(perm.max()),
        "above_both_stratified_permutation_p": float(((perm >= int(both.sum())).sum() + 1) / (perm.size + 1)),
        "above_p1_rr_p1_median": float(np.median(rr_med[1][flag[1]])),
        "above_p1_also_p2_share": float(both.sum() / max(n1, 1)), "above_p2_also_p1_share": float(both.sum() / max(n2, 1)),
        "above_p1_rr_p2_median": float(np.median(rr_med[2][flag[1]])), "above_p1_rr_p2_gt110_share": float((rr_med[2][flag[1]] > 1.10).mean()),
        "above_p1_p2_rr_gt1_share": float((rr_med[2][flag[1]] > 1.0).mean()),
        "above_both_joint_p_ge95": int(((p110[1] * p110[2]) >= 0.95).sum()),
        "above_both_also_pooled": int((both & pooled_flag).sum()),
        "above_either_period": int((flag[1] | flag[2]).sum()), "above_pooled": int(pooled_flag.sum()),
        "above_pooled_above_in_neither_period": int((pooled_flag & ~flag[1] & ~flag[2]).sum()),
        "spearman_rr_p1_p2_all": sp(rr_med[1], rr_med[2]),
        "spearman_rr_p1_p2_ge25_deaths_each": sp(rr_med[1], rr_med[2], big.values), "n_ge25_deaths_each": int(big.sum()),
        "spearman_rr_p1_p2_ge20_expected_each": sp(rr_med[1], rr_med[2], big20.values), "n_ge20_expected_each": int(big20.sum()),
        "spearman_rr_p1_pooled": sp(rr_med[1], pooled_rr), "spearman_rr_p2_pooled": sp(rr_med[2], pooled_rr),
    }
    # regional net excess by period (posterior of sum E_i (RR_i - 1))
    reg = {}
    for p in PERIODS:
        rr = np.exp(res[p]["u"])
        ex = res[p]["E"] * (rr - 1.0)
        reg[p] = {"Brazil": [float(np.median(ex.sum(0))), *map(float, ci(ex.sum(0)))]}
        for r in ["North", "Northeast", "Central-West", "Southeast", "South"]:
            m = (out.REG == r).values
            reg[p][r] = [float(np.median(ex[m].sum(0))), *map(float, ci(ex[m].sum(0)))]
        reg[p]["deaths"] = int(frames[p].deaths_total.sum())
    S["net_excess_by_region"] = reg
    (ROOT / "outputs" / "period_persistence.json").write_text(json.dumps(S, indent=2, default=float))
    print(json.dumps({k: S[k] for k in S if k not in ("net_excess_by_region",)}, indent=1, default=float))


if __name__ == "__main__":
    main()
