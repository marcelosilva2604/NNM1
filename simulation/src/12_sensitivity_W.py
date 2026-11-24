"""Sensitivity of the primary fit (BYM2 + restricted spatial regression) to the neighbourhood matrix.

Step 2 of ANALYSIS_PLAN_2026-11-23 (real data, not simulation). The variants below were chosen when
this script was written (24/11/2026), after the grid finished, and are not in the frozen plan:
    queen_refit  queen contiguity as in the primary fit, refitted here (control: Monte Carlo and
                 code-path noise only, same seed and sampler settings as the primary)
    knn4         4 nearest centroids, symmetrised (union), welded to one component
    knn8         8 nearest centroids, symmetrised (union), welded to one component
Everything else equals the primary fit: spec bym2_rsr_full, 4 chains x (4000 tune + 4000 draws),
target_accept 0.9, seed 20261013. The BYM2 scaling factor is recomputed for each graph.

Comparison with the saved primary posterior (read-only): hyperparameters, and for the primary rule
P(RR > 1.10) >= 0.95 the overlap of flagged municipalities, plus rank agreement of posterior medians.

Run one variant:  python 12_sensitivity_W.py --variant knn4
Compare:          python 12_sensitivity_W.py --compare
"""

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
from scipy import stats

HERE = Path(__file__).resolve().parent
V2 = Path(os.environ.get("NNM1_ROOT", HERE.parent.parent / "primary_fit"))
sys.path.insert(0, str(V2 / "src"))
import nc_io  # noqa: E402

OUT = HERE.parent / "outputs" / "sens_W"
OUT.mkdir(parents=True, exist_ok=True)

_spec = importlib.util.spec_from_file_location("fit_models", V2 / "src" / "04_fit_models.py")
fm = importlib.util.module_from_spec(_spec)
os.environ.setdefault("NNM1_ROOT", str(V2))
_spec.loader.exec_module(fm)
SPEC = fm.SPECS["bym2_rsr_full"]
SEED = 20261013
VARIANTS = ["queen_refit", "knn4", "knn8"]


def components(n, nb):
    lab = np.full(n, -1)
    c = 0
    for s in range(n):
        if lab[s] != -1:
            continue
        st = [s]
        lab[s] = c
        while st:
            a = st.pop()
            for b in nb[a]:
                if lab[b] == -1:
                    lab[b] = c
                    st.append(b)
        c += 1
    return lab


def knn_graph(k):
    from scipy.spatial import cKDTree
    from pyproj import Transformer
    m = pd.read_csv(V2 / "data" / "municipalities.csv")
    x, y = Transformer.from_crs(4326, 5880, always_xy=True).transform(m.lon.values, m.lat.values)
    xy = np.column_stack([x, y])
    n = len(xy)
    _, idx = cKDTree(xy).query(xy, k=k + 1)
    nb = {i: set() for i in range(n)}
    for i in range(n):
        for j in idx[i, 1:]:
            nb[i].add(int(j))
            nb[int(j)].add(i)
    welded = 0
    while True:
        lab = components(n, nb)
        if lab.max() == 0:
            break
        small = int(np.argmin(np.bincount(lab)))
        t = np.where(lab == small)[0]
        o = np.where(lab != small)[0]
        best = None
        for i in t:
            d = np.hypot(xy[o, 0] - xy[i, 0], xy[o, 1] - xy[i, 1])
            q = int(np.argmin(d))
            if best is None or d[q] < best[0]:
                best = (d[q], int(i), int(o[q]))
        _, i, j = best
        nb[i].add(j)
        nb[j].add(i)
        welded += 1
    W = np.zeros((n, n))
    for i, v in nb.items():
        for j in v:
            W[i, j] = 1.0
    return W, welded


def scaling_factor(W):
    Q = np.diag(W.sum(1)) - W
    vals, vecs = np.linalg.eigh(Q)
    keep = vals > 1e-8 * vals.max()
    d = np.einsum("ij,j,ij->i", vecs[:, keep], 1.0 / vals[keep], vecs[:, keep])
    return float(np.exp(np.mean(np.log(d)))), int((~keep).sum())


def run(variant):
    df, Wq, sq = fm.load()
    if variant == "queen_refit":
        W, welded, scaling, zeros = Wq, 0, sq, 1
    else:
        W, welded = knn_graph(int(variant.replace("knn", "")))
        scaling, zeros = scaling_factor(W)
    deg = W.sum(1)
    info = {"variant": variant, "pairs": int(W.sum() / 2), "mean_degree": float(deg.mean()), "min_degree": int(deg.min()),
            "max_degree": int(deg.max()), "welded_joins": welded, "scaling": scaling, "zero_eigenvalues": zeros}
    print(info, flush=True)
    t0 = time.time()
    with fm.build(df, W, scaling, SPEC):
        raw = pm.sample(draws=4000, tune=4000, chains=4, cores=4, random_seed=SEED, target_accept=0.9, progressbar=False)
    p = raw["posterior"].to_dataset()
    rh = az.rhat(p)
    summ = az.summary(raw["posterior"].to_dataset()[["b0", "beta", "sigma_u", "rho"]], ci_prob=0.95, round_to=6)
    info.update({"elapsed_min": round((time.time() - t0) / 60, 1),
                 "divergences": int(raw["sample_stats"].to_dataset()["diverging"].sum()),
                 "rhat_hyper": float(np.nanmax(rh[["b0", "beta", "sigma_u", "rho"]].to_array().values)),
                 "rhat_u": float(np.nanmax(rh["u"].values)),
                 "hyper": summ[["mean", "sd", "eti95_lb", "eti95_ub", "r_hat", "ess_bulk"]].to_dict("index")})
    u = p["u"].stack(s=("chain", "draw")).values
    rr = np.exp(u)
    pd.DataFrame({"CODMUNRES": df.CODMUNRES.values, "births": df.live_births.values, "deaths": df.deaths_total.values,
                  "rr_median": np.median(rr, 1), "p_gt1.10": (rr > 1.10).mean(1), "p_gt1": (rr > 1).mean(1)}
                 ).to_csv(OUT / f"muni_{variant}.csv", index=False)
    (OUT / f"fit_{variant}.json").write_text(json.dumps(info, indent=1, default=str))
    print(f"{variant}: done {info['elapsed_min']} min, div {info['divergences']}, rhat_hyper {info['rhat_hyper']:.3f}, rhat_u {info['rhat_u']:.3f}", flush=True)


def compare():
    df, _, _ = fm.load()
    post = nc_io.load(V2 / "model" / "idata_bym2_rsr_full.nc")["posterior"]
    u = post["u"].stack(s=("chain", "draw")).values
    rr = np.exp(u)
    ref = pd.DataFrame({"CODMUNRES": df.CODMUNRES.values, "rr_median": np.median(rr, 1), "p110": (rr > 1.10).mean(1)})
    ref["flag"] = ref.p110 >= 0.95
    ref_fit = json.loads((V2 / "model" / "fit_bym2_rsr_full.json").read_text())
    rows = [{"variant": "primary (saved)", "flagged": int(ref.flag.sum()), "both": int(ref.flag.sum()), "jaccard": 1.0,
             "spearman_median_rr": 1.0, "max_abs_dp110": 0.0,
             "sigma_u": ref_fit["hyper"]["sigma_u"]["mean"], "rho": ref_fit["hyper"]["rho"]["mean"],
             "beta_nonwhite": ref_fit["hyper"]["beta[sh_nonwhite]"]["mean"], "beta_age_lt20": ref_fit["hyper"]["beta[sh_age_lt20]"]["mean"],
             "divergences": ref_fit["divergences"], "rhat_hyper": ref_fit["max_rhat_hyper"], "rhat_u": ref_fit["max_rhat_u"],
             "mean_degree": 5.80}]
    for v in VARIANTS:
        f = OUT / f"muni_{v}.csv"
        if not f.exists():
            continue
        m = pd.read_csv(f, dtype={"CODMUNRES": str})
        assert (m.CODMUNRES.values == ref.CODMUNRES.values).all()
        fl = (m["p_gt1.10"] >= 0.95).values
        both = int((fl & ref.flag.values).sum())
        info = json.loads((OUT / f"fit_{v}.json").read_text())
        rows.append({"variant": v, "flagged": int(fl.sum()), "both": both,
                     "jaccard": both / max(int((fl | ref.flag.values).sum()), 1),
                     "spearman_median_rr": float(stats.spearmanr(m.rr_median, ref.rr_median)[0]),
                     "max_abs_dp110": float(np.abs(m["p_gt1.10"] - ref.p110).max()),
                     "sigma_u": info["hyper"]["sigma_u"]["mean"], "rho": info["hyper"]["rho"]["mean"],
                     "beta_nonwhite": info["hyper"]["beta[sh_nonwhite]"]["mean"], "beta_age_lt20": info["hyper"]["beta[sh_age_lt20]"]["mean"],
                     "divergences": info["divergences"], "rhat_hyper": info["rhat_hyper"], "rhat_u": info["rhat_u"],
                     "mean_degree": info["mean_degree"]})
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "sensitivity_W_summary.csv", index=False)
    print(d.round(4).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=VARIANTS)
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    compare() if a.compare else run(a.variant)
