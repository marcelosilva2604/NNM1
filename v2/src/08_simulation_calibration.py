"""Simulation-based evaluation of the classification rule (P[RR > 1.10] >= 0.95).

Each replicate takes one posterior draw of the fitted primary model (b0, beta, u) as the known truth,
simulates municipal death counts from Poisson(E_i exp(u_i)) with the observed births and covariates,
refits the same model (BYM2, restricted spatial regression) and compares the classification and the
estimates with the truth. A second scenario simulates with no residual (u = 0) to show how often the
rule flags municipalities when nothing is there. The truth is therefore conditional on a generator
calibrated to the fitted model, and the model is correctly specified: the evaluation shows what the
rule can and cannot do under favourable conditions, not its performance on the real data.

Run one replicate (used with xargs -P to run several in parallel):
    python src/08_simulation_calibration.py --scenario truth --rep 3
Aggregate all replicates written so far:
    python src/08_simulation_calibration.py --aggregate

Needs model/idata_bym2_rsr_full.nc (written by 04). Outputs: outputs/simulation/rep_*.json,
outputs/simulation_summary.json, outputs/tables/simulation_calibration.csv
"""

import argparse
import importlib.util
import json
import os as _os
import time
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

HERE = Path(__file__).resolve().parent
ROOT = Path(_os.environ.get("NNM1_ROOT", HERE.parent))  # repository root; override with NNM1_ROOT
MODEL, DATA = ROOT / "model", ROOT / "data"
SIM = ROOT / "outputs" / "simulation"
TAB = ROOT / "outputs" / "tables"
SIM.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

_spec = importlib.util.spec_from_file_location("fit_models", HERE / "04_fit_models.py")
fm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fm)

SPEC = fm.SPECS["bym2_rsr_full"]
SEED0 = 20261016
STRATA = [("<20", 0, 20), ("20-49", 20, 50), ("50-199", 50, 200), (">=200", 200, np.inf)]
PBINS = [0.0, 0.05, 0.5, 0.8, 0.95, 1.0000001]


def replicate(scenario, k):
    df, W, scaling = fm.load()
    post = az.from_netcdf(MODEL / "idata_bym2_rsr_full.nc")["posterior"]
    u_all = post["u"].stack(s=("chain", "draw")).values
    b0_all = post["b0"].stack(s=("chain", "draw")).values
    beta_all = post["beta"].stack(s=("chain", "draw")).values
    Xr = df[SPEC["X"]].values.astype(float)
    Xz = (Xr - Xr.mean(0)) / Xr.std(0)
    births = df.live_births.values.astype(float)
    rng = np.random.default_rng(SEED0 + (0 if scenario == "truth" else 100000) + k)
    s = int(rng.integers(u_all.shape[1]))
    E_true = births * np.exp(b0_all[s] + Xz @ beta_all[:, s])
    u_true = u_all[:, s] if scenario == "truth" else np.zeros(len(df))
    rr_true = np.exp(u_true)
    y = rng.poisson(E_true * rr_true)
    df2 = df.copy()
    df2["deaths_total"] = y.astype(int)

    t0 = time.time()
    with fm.build(df2, W, scaling, SPEC):
        raw = pm.sample(draws=1000, tune=2000, chains=4, cores=4, random_seed=SEED0 + k + (0 if scenario == "truth" else 5000),
                        target_accept=0.9, progressbar=False)
    p = raw["posterior"].to_dataset()
    u = p["u"].stack(s=("chain", "draw")).values
    b0 = p["b0"].stack(s=("chain", "draw")).values
    beta = p["beta"].stack(s=("chain", "draw")).values
    rr = np.exp(u)
    p110 = (rr > 1.10).mean(1)
    flag = p110 >= 0.95
    lo, hi = np.quantile(rr, [0.025, 0.975], axis=1)
    gt110, gt125 = rr_true > 1.10, rr_true > 1.25

    def strata_stats(mask_true, name):
        out = {}
        for lab, a, b in STRATA:
            m = (E_true >= a) & (E_true < b)
            out[lab] = {"n_true": int((mask_true & m).sum()), "n_flagged_of_true": int((mask_true & m & flag).sum())}
        return out

    cal = []
    for a, b in zip(PBINS[:-1], PBINS[1:]):
        m = (p110 >= a) & (p110 < b)
        cal.append({"lo": a, "hi": min(b, 1.0), "n": int(m.sum()), "n_true_gt110": int((m & gt110).sum())})
    cover = (lo <= rr_true) & (rr_true <= hi)
    cover_s = {lab: {"n": int(((E_true >= a) & (E_true < b)).sum()), "n_covered": int((cover & (E_true >= a) & (E_true < b)).sum())}
               for lab, a, b in STRATA}
    E_hat_all = births[:, None] * np.exp(b0[None, :] + Xz @ beta)
    ex_all = (E_hat_all * (rr - 1.0)).sum(0)
    if flag.any():
        ex_flag = (E_hat_all[flag] * (rr[flag] - 1.0)).sum(0)
        est_flag = float(np.median(ex_flag))
        est_flag_lo, est_flag_hi = map(float, np.quantile(ex_flag, [0.025, 0.975]))
    else:
        est_flag = est_flag_lo = est_flag_hi = float("nan")
    true_flag = float((E_true[flag] * (rr_true[flag] - 1.0)).sum())
    res = {
        "scenario": scenario, "rep": k, "draw": s, "elapsed_min": round((time.time() - t0) / 60, 1),
        "divergences": int(raw["sample_stats"].to_dataset()["diverging"].sum()),
        "n_flag": int(flag.sum()), "tp": int((flag & gt110).sum()), "fp": int((flag & ~gt110).sum()),
        "n_true_gt110": int(gt110.sum()), "n_true_gt125": int(gt125.sum()),
        "flagged_of_true_gt110": int((flag & gt110).sum()), "flagged_of_true_gt125": int((flag & gt125).sum()),
        "strata_gt110": strata_stats(gt110, "gt110"), "strata_gt125": strata_stats(gt125, "gt125"),
        "coverage_overall": float(cover.mean()), "coverage_strata": cover_s, "calibration": cal,
        "flagged_below": int((((rr < 0.90).mean(1)) >= 0.95).sum()),
        "true_lt090": int((rr_true < 0.90).sum()),
        "net_excess_true": float((E_true * (rr_true - 1.0)).sum()),
        "net_excess_est_median": float(np.median(ex_all)),
        "net_excess_est_lo95": float(np.quantile(ex_all, 0.025)), "net_excess_est_hi95": float(np.quantile(ex_all, 0.975)),
        "flagged_excess_true": true_flag, "flagged_excess_est_median": est_flag,
        "flagged_excess_est_lo95": est_flag_lo, "flagged_excess_est_hi95": est_flag_hi,
        "sd_u_true": float(np.std(u_true)), "sd_u_est_median": float(np.median(u.std(0))),
    }
    (SIM / f"rep_{scenario}_{k}.json").write_text(json.dumps(res, indent=1))
    print(f"{scenario} {k}: flagged {res['n_flag']}, fp {res['fp']}, true>1.10 {res['n_true_gt110']}, {res['elapsed_min']} min", flush=True)


def aggregate():
    reps = [json.loads(p.read_text()) for p in sorted(SIM.glob("rep_*.json"))]
    T = [r for r in reps if r["scenario"] == "truth"]
    N = [r for r in reps if r["scenario"] == "null"]
    S = {"n_truth_replicates": len(T), "n_null_replicates": len(N)}
    if T:
        flagged = np.array([r["n_flag"] for r in T]); fp = np.array([r["fp"] for r in T])
        S["truth"] = {
            "flagged_mean": float(flagged.mean()), "flagged_min": int(flagged.min()), "flagged_max": int(flagged.max()),
            "true_gt110_mean": float(np.mean([r["n_true_gt110"] for r in T])),
            "false_flag_share_pooled": float(fp.sum() / flagged.sum()),
            "false_flag_share_by_rep_min": float(np.min(fp / np.maximum(flagged, 1))),
            "false_flag_share_by_rep_max": float(np.max(fp / np.maximum(flagged, 1))),
            "sensitivity_gt110_pooled": float(sum(r["flagged_of_true_gt110"] for r in T) / sum(r["n_true_gt110"] for r in T)),
            "sensitivity_gt125_pooled": float(sum(r["flagged_of_true_gt125"] for r in T) / sum(r["n_true_gt125"] for r in T)),
            "coverage_overall_mean": float(np.mean([r["coverage_overall"] for r in T])),
            "net_excess_bias_mean": float(np.mean([r["net_excess_est_median"] - r["net_excess_true"] for r in T])),
            "net_excess_true_mean": float(np.mean([r["net_excess_true"] for r in T])),
            "net_excess_ci_covers_truth": float(np.mean([r["net_excess_est_lo95"] <= r["net_excess_true"] <= r["net_excess_est_hi95"] for r in T])),
            "flagged_excess_overstatement_mean": float(np.nanmean([r["flagged_excess_est_median"] - r["flagged_excess_true"] for r in T])),
            "flagged_excess_true_mean": float(np.nanmean([r["flagged_excess_true"] for r in T])),
            "flagged_excess_est_mean": float(np.nanmean([r["flagged_excess_est_median"] for r in T])),
            "flagged_excess_ci_covers_truth": float(np.nanmean([r["flagged_excess_est_lo95"] <= r["flagged_excess_true"] <= r["flagged_excess_est_hi95"] for r in T])),
            "divergences_total": int(sum(r["divergences"] for r in T)),
        }
        rows = []
        for lab, _, _ in STRATA:
            t110 = sum(r["strata_gt110"][lab]["n_true"] for r in T); f110 = sum(r["strata_gt110"][lab]["n_flagged_of_true"] for r in T)
            t125 = sum(r["strata_gt125"][lab]["n_true"] for r in T); f125 = sum(r["strata_gt125"][lab]["n_flagged_of_true"] for r in T)
            nn = sum(r["coverage_strata"][lab]["n"] for r in T); cc = sum(r["coverage_strata"][lab]["n_covered"] for r in T)
            rows.append({"expected_deaths": lab, "n_true_gt110": t110, "flagged_of_true_gt110": f110,
                         "sensitivity_gt110": f110 / t110 if t110 else np.nan, "n_true_gt125": t125, "flagged_of_true_gt125": f125,
                         "sensitivity_gt125": f125 / t125 if t125 else np.nan, "coverage_95": cc / nn if nn else np.nan})
        pd.DataFrame(rows).to_csv(TAB / "simulation_calibration.csv", index=False)
        cal = []
        for i, (a, b) in enumerate(zip(PBINS[:-1], PBINS[1:])):
            n_ = sum(r["calibration"][i]["n"] for r in T); t_ = sum(r["calibration"][i]["n_true_gt110"] for r in T)
            cal.append({"p_lo": a, "p_hi": min(b, 1.0), "n": n_, "share_truly_gt110": t_ / n_ if n_ else np.nan})
        pd.DataFrame(cal).to_csv(TAB / "simulation_calibration_probability.csv", index=False)
        S["strata"] = rows; S["probability_calibration"] = cal
    if N:
        nf = np.array([r["n_flag"] for r in N]); nb = np.array([r["flagged_below"] for r in N])
        S["null"] = {"flagged_above_mean": float(nf.mean()), "flagged_above_max": int(nf.max()), "flagged_below_mean": float(nb.mean()),
                     "reps_with_any_false_flag": int((nf > 0).sum())}
    (ROOT / "outputs" / "simulation_summary.json").write_text(json.dumps(S, indent=2, default=float))
    print(json.dumps({k: v for k, v in S.items() if k not in ("strata", "probability_calibration")}, indent=1, default=float))
    if "strata" in S:
        print(pd.DataFrame(S["strata"]).round(3).to_string(index=False))
        print(pd.DataFrame(S["probability_calibration"]).round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", choices=["truth", "null"])
    ap.add_argument("--rep", type=int)
    ap.add_argument("--aggregate", action="store_true")
    a = ap.parse_args()
    if a.aggregate:
        aggregate()
    else:
        replicate(a.scenario, a.rep)
