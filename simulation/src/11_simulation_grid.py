"""Grid simulation for the SMMR version: operating characteristics of posterior-probability rules.

One replicate = one posterior draw of the fitted primary model (b0, beta, u) as known truth, with the
municipal effect scaled by k, Poisson counts at the observed births, a full refit of the primary model
(BYM2 + restricted spatial regression), then evaluation of a family of rules and three comparators from
the same posterior (no further fitting). Plan frozen in ../ANALYSIS_PLAN_2026-11-23.md.

Run one replicate (used with xargs -P):
    python 11_simulation_grid.py --k 1.0 --rep 3
Aggregate:
    python 11_simulation_grid.py --aggregate

Uses the primary fit in ../../primary_fit (model/idata_bym2_rsr_full.nc, data/*.csv), read-only; set NNM1_ROOT to override.
"""

import argparse
import importlib.util
import json
import os
import time
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
from scipy import stats

HERE = Path(__file__).resolve().parent
# Primary fit (the posterior used as generating truth). Read-only: never write inside this folder.
V2 = Path(os.environ.get("NNM1_ROOT", HERE.parent.parent / "primary_fit"))
import sys
sys.path.insert(0, str(V2 / "src"))
import nc_io
OUT = HERE.parent / "outputs" / "grid"
OUT.mkdir(parents=True, exist_ok=True)

_spec = importlib.util.spec_from_file_location("fit_models", V2 / "src" / "04_fit_models.py")
fm = importlib.util.module_from_spec(_spec)
os.environ.setdefault("NNM1_ROOT", str(V2))
_spec.loader.exec_module(fm)

SPEC = fm.SPECS["bym2_rsr_full"]
SEED0 = 20261123
KS = [0.0, 0.5, 1.0, 1.5, 2.0]
NULL_REPS = 40  # amendment 24/11/2026: null scenario limited to 40 replicates, extended sampling (see ANALYSIS_PLAN amendment)
TAUS = [1.05, 1.10, 1.20]
CUTS = [0.80, 0.90, 0.95, 0.99]
STRATA = [("<20", 0, 20), ("20-49", 20, 50), ("50-199", 50, 200), (">=200", 200, np.inf)]
PBINS = [0.0, 0.05, 0.5, 0.8, 0.95, 1.0000001]


def bh_flag(p, q=0.05):
    n = len(p)
    order = np.argsort(p)
    thresh = q * (np.arange(1, n + 1) / n)
    ok = p[order] <= thresh
    flag = np.zeros(n, bool)
    if ok.any():
        flag[order[: np.max(np.where(ok)[0]) + 1]] = True
    return flag


def evaluate(flag, truth, E_true):
    d = {"n_flag": int(flag.sum()), "tp": int((flag & truth).sum()), "fp": int((flag & ~truth).sum()),
         "n_true": int(truth.sum())}
    d["strata"] = {}
    for lab, a, b in STRATA:
        m = (E_true >= a) & (E_true < b)
        d["strata"][lab] = {"n_true": int((truth & m).sum()), "tp": int((flag & truth & m).sum()),
                            "n_flag": int((flag & m).sum()), "fp": int((flag & ~truth & m).sum())}
    return d


def replicate(k, rep):
    if k == 0.0 and rep >= NULL_REPS:
        return
    draws, tune = (4000, 4000) if k == 0.0 else (1000, 2000)
    ki = KS.index(k)
    df, W, scaling = fm.load()
    post = nc_io.load(V2 / "model" / "idata_bym2_rsr_full.nc")["posterior"]
    u_all = post["u"].stack(s=("chain", "draw")).values
    b0_all = post["b0"].stack(s=("chain", "draw")).values
    beta_all = post["beta"].stack(s=("chain", "draw")).values
    Xr = df[SPEC["X"]].values.astype(float)
    Xz = (Xr - Xr.mean(0)) / Xr.std(0)
    births = df.live_births.values.astype(float)
    seed = SEED0 + 1000 * ki + rep
    rng = np.random.default_rng(seed)
    s = int(rng.integers(u_all.shape[1]))
    E_true = births * np.exp(b0_all[s] + Xz @ beta_all[:, s])
    u_true = k * u_all[:, s]
    rr_true = np.exp(u_true)
    y = rng.poisson(E_true * rr_true)
    df2 = df.copy()
    df2["deaths_total"] = y.astype(int)

    t0 = time.time()
    with fm.build(df2, W, scaling, SPEC):
        raw = pm.sample(draws=draws, tune=tune, chains=4, cores=4, random_seed=seed, target_accept=0.9, progressbar=False)
    p = raw["posterior"].to_dataset()
    u = p["u"].stack(s=("chain", "draw")).values
    b0 = p["b0"].stack(s=("chain", "draw")).values
    beta = p["beta"].stack(s=("chain", "draw")).values
    rr = np.exp(u)
    lo, hi = np.quantile(rr, [0.025, 0.975], axis=1)
    cover = (lo <= rr_true) & (rr_true <= hi)

    rules = {}
    for tau in TAUS:
        pt = (rr > tau).mean(1)
        truth = rr_true > tau
        for c in CUTS:
            rules[f"post_tau{tau:.2f}_c{c:.2f}"] = evaluate(pt >= c, truth, E_true)
    truth110 = rr_true > 1.10
    rules["C1_median_gt1.10"] = evaluate(np.median(rr, axis=1) > 1.10, truth110, E_true)
    E_hat = (births[:, None] * np.exp(b0[None, :] + Xz @ beta)).mean(1)
    pval = stats.poisson.sf(y - 1, 1.10 * E_hat)  # P(Y >= y | mean = 1.10 E_hat)
    rules["C2_poisson_p05"] = evaluate(pval < 0.05, truth110, E_true)
    rules["C3_poisson_BH05"] = evaluate(bh_flag(pval, 0.05), truth110, E_true)

    p110 = (rr > 1.10).mean(1)
    flag_pri = p110 >= 0.95
    cal = [{"lo": a, "hi": min(b, 1.0), "n": int(((p110 >= a) & (p110 < b)).sum()),
            "n_true_gt110": int((((p110 >= a) & (p110 < b)) & truth110).sum())} for a, b in zip(PBINS[:-1], PBINS[1:])]
    cov_strata = {lab: {"n": int(((E_true >= a) & (E_true < b)).sum()),
                        "n_cov": int((cover & (E_true >= a) & (E_true < b)).sum())} for lab, a, b in STRATA}
    rh = az.rhat(raw["posterior"].to_dataset())
    rhat = float(np.nanmax(rh.to_array().values))
    rhat_hyper = float(np.nanmax(rh[["b0", "beta", "sigma_u", "rho"]].to_array().values))
    rhat_u = float(np.nanmax(rh["u"].values))
    res = {"k": k, "rep": rep, "seed": seed, "draw": s, "elapsed_min": round((time.time() - t0) / 60, 1),
           "divergences": int(raw["sample_stats"].to_dataset()["diverging"].sum()), "rhat_max": rhat, "rhat_hyper": rhat_hyper, "rhat_u": rhat_u, "draws": draws, "tune": tune,
           "sd_u_true": float(np.std(u_true)), "sd_u_est_median": float(np.median(u.std(0))),
           "n_true_gt110": int(truth110.sum()), "coverage_overall": float(cover.mean()),
           "coverage_flagged": float(cover[flag_pri].mean()) if flag_pri.any() else None,
           "coverage_strata": cov_strata, "calibration": cal, "rules": rules}
    (OUT / f"rep_k{k:.1f}_{rep}.json").write_text(json.dumps(res))
    print(f"k={k} rep={rep}: flagged(primary) {rules['post_tau1.10_c0.95']['n_flag']}, fp {rules['post_tau1.10_c0.95']['fp']}, "
          f"true>1.10 {int(truth110.sum())}, {res['elapsed_min']} min, div {res['divergences']}, rhat {rhat:.3f}", flush=True)


def aggregate():
    reps = [json.loads(p.read_text()) for p in sorted(OUT.glob("rep_k*_*.json"))]
    rows = []
    for k in KS:
        R = [r for r in reps if r["k"] == k]
        if not R:
            continue
        for rule in R[0]["rules"]:
            for lab in ["ALL"] + [s[0] for s in STRATA]:
                def g(r, key):
                    d = r["rules"][rule]
                    return d[key] if lab == "ALL" else d["strata"][lab][key]
                nt = np.array([g(r, "n_true") for r in R], float)
                tp = np.array([g(r, "tp") for r in R], float)
                nf = np.array([g(r, "n_flag") for r in R], float)
                fp = np.array([g(r, "fp") for r in R], float)
                sens_rep = np.where(nt > 0, tp / np.maximum(nt, 1), np.nan)
                fdp_rep = np.where(nf > 0, fp / np.maximum(nf, 1), np.nan)
                rows.append({"k": k, "rule": rule, "stratum": lab, "n_reps": len(R),
                             "n_true_total": int(nt.sum()), "n_flag_total": int(nf.sum()),
                             "sens_pooled": tp.sum() / nt.sum() if nt.sum() else np.nan,
                             "sens_rep_mean": np.nanmean(sens_rep) if np.isfinite(sens_rep).any() else np.nan,
                             "sens_rep_mcse": np.nanstd(sens_rep, ddof=1) / np.sqrt(np.isfinite(sens_rep).sum()) if np.isfinite(sens_rep).sum() > 1 else np.nan,
                             "false_flag_pooled": fp.sum() / nf.sum() if nf.sum() else np.nan,
                             "false_flag_rep_mean": np.nanmean(fdp_rep) if np.isfinite(fdp_rep).any() else np.nan,
                             "false_flag_rep_mcse": np.nanstd(fdp_rep, ddof=1) / np.sqrt(np.isfinite(fdp_rep).sum()) if np.isfinite(fdp_rep).sum() > 1 else np.nan,
                             "reps_with_any_flag": int((nf > 0).sum())})
    pd.DataFrame(rows).to_csv(OUT.parent / "grid_summary.csv", index=False)
    meta = [{"k": r["k"], "rep": r["rep"], "elapsed_min": r["elapsed_min"], "divergences": r["divergences"], "rhat_max": r["rhat_max"], "rhat_hyper": r.get("rhat_hyper"), "rhat_u": r.get("rhat_u"), "draws": r.get("draws", 1000),
             "sd_u_true": r["sd_u_true"], "sd_u_est_median": r["sd_u_est_median"], "coverage_overall": r["coverage_overall"],
             "coverage_flagged": r["coverage_flagged"]} for r in reps]
    pd.DataFrame(meta).to_csv(OUT.parent / "grid_replicate_meta.csv", index=False)
    cal = []
    for k in KS:
        R = [r for r in reps if r["k"] == k]
        for i in range(len(PBINS) - 1):
            n_ = sum(r["calibration"][i]["n"] for r in R)
            t_ = sum(r["calibration"][i]["n_true_gt110"] for r in R)
            if R:
                cal.append({"k": k, "p_lo": PBINS[i], "p_hi": min(PBINS[i + 1], 1.0), "n": n_, "share_truly_gt110": t_ / n_ if n_ else np.nan})
    pd.DataFrame(cal).to_csv(OUT.parent / "grid_calibration.csv", index=False)
    print(len(reps), "replicates aggregated")
    d = pd.DataFrame(rows)
    print(d[(d.rule == "post_tau1.10_c0.95") & (d.stratum == "ALL")].round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=float)
    ap.add_argument("--rep", type=int)
    ap.add_argument("--aggregate", action="store_true")
    a = ap.parse_args()
    if a.aggregate:
        aggregate()
    else:
        replicate(a.k, a.rep)
