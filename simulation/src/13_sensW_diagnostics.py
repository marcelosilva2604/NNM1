"""Diagnostics for the neighbourhood-matrix sensitivity (12_sensitivity_W.py). Investigation only.

Two questions:
  1. Why does sigma_u fall from about 0.37 (queen) to about 0.19 (kNN graphs)?
  2. Where are the divergences of the knn8 fit, and do they go away with a higher target_accept?

Nothing written by 12_sensitivity_W.py is modified. The graphs are built by importing that script, so
the objects checked here are the ones it used. Outputs go to outputs/sens_W/diag/.

Modes:
    --graphs                          graph checks, no sampling (order, symmetry, degree, kNN rebuilt by
                                      brute force, BYM2 scaling recomputed two ways, spectrum of Q)
    --fit VARIANT --target-accept X   refit one variant (same model, seed, 4 chains x 4000/4000) keeping
                                      phi and theta in memory to locate divergences and to measure the
                                      realised spread of u and of its spatial part
    --report                          gather the diag JSON files into one table

Run through ./run.sh, for example:
    ./run.sh src/13_sensW_diagnostics.py --graphs
    ./run.sh src/13_sensW_diagnostics.py --fit knn8 --target-accept 0.97
"""

import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DIAG = HERE.parent / "outputs" / "sens_W" / "diag"
DIAG.mkdir(parents=True, exist_ok=True)
SENS = HERE.parent / "outputs" / "sens_W"

_spec = importlib.util.spec_from_file_location("sens12", HERE / "12_sensitivity_W.py")
s12 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s12)
fm, V2, SEED, SPEC = s12.fm, s12.V2, s12.SEED, s12.SPEC
VARIANTS = ["queen", "knn4", "knn8"]
BANDS = [(0, 10), (10, 100), (100, 1000), (1000, None)]  # eigen-index bands, smoothest modes first


def get_graph(variant):
    """Adjacency and scaling exactly as 12_sensitivity_W.py obtains them."""
    df, Wq, sq = fm.load()
    if variant == "queen":
        return df, Wq, float(sq), 0
    W, welded = s12.knn_graph(int(variant.replace("knn", "")))
    scaling, _ = s12.scaling_factor(W)
    return df, W, float(scaling), welded


def spectrum(W):
    Q = np.diag(W.sum(1)) - W
    vals, vecs = np.linalg.eigh(Q)
    keep = vals > 1e-8 * vals.max()
    return vals, vecs, keep


def projector(df):
    X = df[SPEC["X"]].values.astype(float)
    Xz = (X - X.mean(0)) / X.std(0)
    Xa = np.column_stack([np.ones(len(df)), Xz])
    return np.eye(len(df)) - Xa @ np.linalg.inv(Xa.T @ Xa) @ Xa.T


def brute_knn(xy, k):
    """Symmetrised (union) kNN adjacency by brute-force distances, independent of the KD-tree."""
    n = len(xy)
    W = np.zeros((n, n))
    ties = 0
    for a in range(0, n, 500):
        b = min(a + 500, n)
        d = np.hypot(xy[a:b, None, 0] - xy[None, :, 0], xy[a:b, None, 1] - xy[None, :, 1])
        d[np.arange(b - a), np.arange(a, b)] = np.inf
        part = np.argpartition(d, k, axis=1)[:, :k + 1]
        ds = np.sort(np.take_along_axis(d, part, 1), axis=1)
        ties += int((ds[:, k] == ds[:, k - 1]).sum())
        idx = np.argsort(d, axis=1)[:, :k]
        for r in range(b - a):
            W[a + r, idx[r]] = 1.0
    return np.maximum(W, W.T), ties


def graphs():
    from pyproj import Transformer
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    df = pd.read_csv(V2 / "data" / "analysis_muni.csv", dtype={"CODMUNRES": str})
    m = pd.read_csv(V2 / "data" / "municipalities.csv", dtype={"CODMUNRES": str})
    order_ok = bool(len(m) == len(df) and (m.CODMUNRES.values == df.CODMUNRES.values).all())
    x, y = Transformer.from_crs(4326, 5880, always_xy=True).transform(m.lon.values, m.lat.values)
    xy = np.column_stack([x, y])
    dup_points = int(len(xy) - len(np.unique(xy, axis=0)))
    stored_queen = json.loads((V2 / "data" / "bym2_scaling.json").read_text())
    rows, per_muni = [], {"CODMUNRES": df.CODMUNRES.values}
    for v in VARIANTS:
        t0 = time.time()
        dfv, W, scaling_used, welded = get_graph(v)
        n = len(W)
        P = projector(dfv)
        deg = W.sum(1)
        ncomp, _ = connected_components(csr_matrix(W), directed=False)
        vals, vecs, keep = spectrum(W)
        inv = 1.0 / vals[keep]
        d = np.einsum("ij,j,ij->i", vecs[:, keep], inv, vecs[:, keep])
        scaling_02 = float(np.exp(np.mean(np.log(d))))
        # Second route to diag(Q^+): (Q + 11'/n)^-1 - 11'/n, no eigen-decomposition.
        Q = np.diag(deg) - W
        d2 = np.diag(np.linalg.inv(Q + 1.0 / n)) - 1.0 / n
        scaling_inv = float(np.exp(np.mean(np.log(d2))))
        # Prior marginal variance of the restricted (projected) scaled ICAR component, sigma_u = rho = 1.
        A = P @ vecs[:, keep]
        d_rsr = np.einsum("ij,j,ij->i", A, inv, A) / scaling_02
        tr = inv.sum()
        row = {
            "variant": v, "n": n, "order_municipalities_matches_analysis_table": order_ok,
            "symmetric": bool((W == W.T).all()), "binary": bool(np.isin(W, [0.0, 1.0]).all()),
            "zero_diagonal": bool((np.diag(W) == 0).all()), "components": int(ncomp), "welded_joins": welded,
            "pairs": int(W.sum() / 2), "mean_degree": float(deg.mean()), "min_degree": int(deg.min()),
            "max_degree": int(deg.max()), "zero_eigenvalues": int((~keep).sum()),
            "scaling_used_by_script12": scaling_used, "scaling_formula_02": scaling_02, "scaling_via_inverse": scaling_inv,
            "rel_diff_used_vs_02": abs(scaling_used - scaling_02) / scaling_02,
            "diagQplus_min": float(d.min()), "diagQplus_median": float(np.median(d)), "diagQplus_max": float(d.max()),
            "diagQplus_arith_mean": float(d.mean()), "arith_over_geo_mean": float(d.mean() / scaling_02),
            "lambda2": float(vals[keep][0]), "lambda_max": float(vals.max()),
            "prior_var_rsr_spatial_geo_mean": float(np.exp(np.mean(np.log(d_rsr)))),
            "prior_var_rsr_spatial_arith_mean": float(d_rsr.mean()),
        }
        for a, b in BANDS:
            row[f"share_trace_Qplus_modes_{a}_{b or 'end'}"] = float(inv[a:b].sum() / tr)
        if v == "queen":
            row["scaling_stored_json"] = stored_queen["scaling_factor"]
            row["pairs_stored_json"] = stored_queen["pairs"]
        else:
            k = int(v.replace("knn", ""))
            Wb, ties = brute_knn(xy, k)
            row["knn_equals_bruteforce"] = bool((Wb == W).all())
            row["knn_edges_differing_from_bruteforce"] = int(np.abs(Wb - W).sum() / 2)
            row["knn_distance_ties_at_k"] = ties
            row["duplicate_centroids"] = dup_points
            used = json.loads((SENS / f"fit_{v}.json").read_text())
            row["scaling_in_fit_json"] = used["scaling"]
            row["pairs_in_fit_json"] = used["pairs"]
            # Mean neighbour distance, to describe the spatial reach of each graph.
        i, j = np.where(np.triu(W) == 1)
        dist = np.hypot(xy[i, 0] - xy[j, 0], xy[i, 1] - xy[j, 1]) / 1000.0
        row["edge_km_median"] = float(np.median(dist))
        row["edge_km_p99"] = float(np.quantile(dist, 0.99))
        row["edge_km_max"] = float(dist.max())
        row["elapsed_s"] = round(time.time() - t0, 1)
        rows.append(row)
        per_muni[f"degree_{v}"] = deg.astype(int)
        per_muni[f"diagQplus_{v}"] = d
        print(json.dumps(row, indent=1), flush=True)
    # Edge overlap between graphs.
    Ws = {v: get_graph(v)[1] for v in VARIANTS}
    for a in VARIANTS:
        for b in VARIANTS:
            if a < b:
                both = int((Ws[a] * Ws[b]).sum() / 2)
                print(f"edges shared {a} & {b}: {both} (of {int(Ws[a].sum()/2)} and {int(Ws[b].sum()/2)})", flush=True)
                rows.append({"variant": f"overlap_{a}_{b}", "pairs": both})
    pd.DataFrame(rows).to_csv(DIAG / "graph_checks.csv", index=False)
    (DIAG / "graph_checks.json").write_text(json.dumps(rows, indent=1))
    pd.DataFrame(per_muni).to_csv(DIAG / "graph_per_muni.csv", index=False)


def q3(a):
    return {"mean": float(np.mean(a)), "q025": float(np.quantile(a, 0.025)), "q50": float(np.quantile(a, 0.5)),
            "q975": float(np.quantile(a, 0.975))}


def fit(variant, target_accept, suffix=""):
    import arviz as az
    import pymc as pm

    tag = f"{variant}_ta{target_accept:g}{suffix}"
    df, W, scaling, welded = get_graph(variant)
    n = len(df)
    deg = W.sum(1)
    info = {"variant": variant, "target_accept": target_accept, "seed": SEED, "chains": 4, "tune": 4000, "draws": 4000,
            "scaling": scaling, "pairs": int(W.sum() / 2), "mean_degree": float(deg.mean())}
    print(info, flush=True)
    t0 = time.time()
    with fm.build(df, W, scaling, SPEC):
        raw = pm.sample(draws=4000, tune=4000, chains=4, cores=4, random_seed=SEED, target_accept=target_accept,
                        progressbar=False)
    info["elapsed_min"] = round((time.time() - t0) / 60, 1)
    p = raw["posterior"].to_dataset()
    ss = raw["sample_stats"].to_dataset()
    hyp = p[["b0", "beta", "sigma_u", "rho"]]
    rh = az.rhat(p[["b0", "beta", "sigma_u", "rho", "u"]])
    summ = az.summary(hyp, ci_prob=0.95, round_to=6)
    div = ss["diverging"].values.astype(bool)  # chain x draw
    info.update({
        "divergences": int(div.sum()), "divergences_per_chain": [int(c) for c in div.sum(1)],
        "rhat_hyper": float(np.nanmax(rh[["b0", "beta", "sigma_u", "rho"]].to_array().values)),
        "rhat_u": float(np.nanmax(rh["u"].values)),
        "ess_bulk_min_u": float(az.ess(p[["u"]])["u"].min()),
        "ess_tail_sigma_u": float(az.ess(p[["sigma_u"]], method="tail")["sigma_u"]),
        "ess_tail_rho": float(az.ess(p[["rho"]], method="tail")["rho"]),
        "hyper": summ[["mean", "sd", "eti95_lb", "eti95_ub", "r_hat", "ess_bulk"]].to_dict("index"),
    })
    for name in ["step_size", "tree_depth", "n_steps", "reached_max_treedepth"]:
        if name in ss:
            info[f"{name}_mean_per_chain"] = [float(c) for c in ss[name].values.astype(float).mean(1)]

    flat = lambda a: a.reshape((-1,) + a.shape[2:])
    sig, rho = flat(p["sigma_u"].values), flat(p["rho"].values)
    phi, theta, u = flat(p["phi"].values), flat(p["theta"].values), flat(p["u"].values)
    d = flat(div)
    S = len(sig)

    # Realised spread across municipalities, per draw. These are on the same scale for every graph,
    # unlike sigma_u, whose meaning depends on the graph through the scaling factor.
    spatial = (sig * np.sqrt(rho / scaling))[:, None] * phi
    iid = (sig * np.sqrt(1.0 - rho))[:, None] * theta
    edges_i, edges_j = np.where(np.tril(W) == 1)
    pqp = np.empty(S)
    for a in range(0, S, 2000):
        pqp[a:a + 2000] = ((phi[a:a + 2000][:, edges_i] - phi[a:a + 2000][:, edges_j]) ** 2).sum(1)
    per_draw = {
        "sigma_u": sig, "rho": rho, "sd_u_across_muni": u.std(1), "sd_spatial_across_muni": spatial.std(1),
        "sd_iid_across_muni": iid.std(1), "tau_inv_sigma2rho_over_scaling": sig ** 2 * rho / scaling,
        "iid_var_sigma2_1mrho": sig ** 2 * (1.0 - rho), "phiQphi_over_nm1": pqp / (n - 1),
        "edge_diff_sd_spatial": np.sqrt(sig ** 2 * rho / scaling * pqp / len(edges_i)),
    }
    info["realised"] = {k: q3(v) for k, v in per_draw.items()}
    u_mean = u.mean(0)
    info["sd_of_posterior_mean_u"] = float(u_mean.std())
    info["sd_of_posterior_mean_spatial"] = float(spatial.mean(0).std())
    info["mean_posterior_sd_u"] = float(u.std(0).mean())

    # Spectrum: second moment of phi along each eigenvector of Q, times the eigenvalue. Under the ICAR
    # prior with sigma = 1 this is 1 for every mode; the posterior value shows which spatial scales
    # carry less (or more) variance than the prior implies.
    vals, vecs, keep = spectrum(W)
    lam, Vk = vals[keep], vecs[:, keep]
    m2 = np.zeros(Vk.shape[1])
    for a in range(0, S, 1000):
        m2 += ((phi[a:a + 1000] @ Vk) ** 2).sum(0)
    m2 /= S
    ratio = m2 * lam
    tau_inv = float(per_draw["tau_inv_sigma2rho_over_scaling"].mean())
    bands = []
    for a, b in BANDS:
        sl = slice(a, b)
        bands.append({"modes": f"{a}-{b or len(lam)}", "prior_share_of_trace_Qplus": float((1 / lam[sl]).sum() / (1 / lam).sum()),
                      "posterior_over_prior_second_moment_mean": float(ratio[sl].mean()),
                      "posterior_over_prior_second_moment_median": float(np.median(ratio[sl])),
                      "realised_mean_variance_contribution": float(tau_inv * m2[sl].sum() / n),
                      "prior_mean_variance_contribution_at_posterior_tau": float(tau_inv * (1 / lam[sl]).sum() / n)})
    info["spectrum_bands"] = bands
    info["realised_spatial_mean_variance_all_modes"] = float(tau_inv * m2.sum() / n)
    info["prior_spatial_mean_variance_all_modes_at_posterior_tau"] = float(tau_inv * (1 / lam).sum() / n)

    # Divergences: where in hyperparameter space, and which municipalities.
    muni = pd.DataFrame({"CODMUNRES": df.CODMUNRES.values, "births": df.live_births.values, "deaths": df.deaths_total.values,
                         "degree": deg.astype(int), "u_mean": u_mean, "u_sd": u.std(0), "rr_median": np.median(np.exp(u), 0),
                         "p_gt1.10": (np.exp(u) > 1.10).mean(0), "spatial_mean": spatial.mean(0), "iid_mean": iid.mean(0),
                         "phi_sd": phi.std(0), "theta_sd": theta.std(0)})
    if d.sum() > 0:
        dv = {}
        for name, a in [("rho", rho), ("sigma_u", sig), ("sqrt_1_minus_rho", np.sqrt(1 - rho)),
                        ("sd_iid_across_muni", per_draw["sd_iid_across_muni"]),
                        ("theta_norm_over_sqrt_n", np.sqrt((theta ** 2).mean(1)))]:
            dv[name] = {"divergent": q3(a[d]), "non_divergent": q3(a[~d])}
        for cut in [0.95, 0.98, 0.99]:
            dv[f"share_draws_rho_gt_{cut}"] = {"divergent": float((rho[d] > cut).mean()), "non_divergent": float((rho[~d] > cut).mean())}
        edges = np.quantile(rho, np.linspace(0, 1, 11))
        dv["divergence_rate_by_rho_decile"] = [
            {"rho_lo": float(edges[k]), "rho_hi": float(edges[k + 1]),
             "rate": float(d[(rho >= edges[k]) & ((rho < edges[k + 1]) if k < 9 else (rho <= edges[k + 1]))].mean())}
            for k in range(10)]
        for name in ["energy", "lp", "tree_depth", "step_size"]:
            if name in ss:
                a = flat(ss[name].values.astype(float))
                dv[name] = {"divergent": q3(a[d]), "non_divergent": q3(a[~d])}
        # Standardised mean difference, divergent minus non-divergent draws, per municipality.
        for name, a in [("phi", phi), ("theta", theta), ("u", u)]:
            z = (a[d].mean(0) - a[~d].mean(0)) / a.std(0)
            muni[f"smd_div_{name}"] = z
            dv[f"smd_{name}_abs_max"] = float(np.abs(z).max())
            dv[f"smd_{name}_abs_q99"] = float(np.quantile(np.abs(z), 0.99))
            dv[f"smd_{name}_abs_median"] = float(np.median(np.abs(z)))
        # Reference for the SMD: the same statistic with the divergent labels shifted by half the run
        # within each chain (keeps the count and the clustering of the labels, breaks the link).
        dshift = flat(np.roll(div, div.shape[1] // 2, axis=1))
        for name, a in [("phi", phi), ("theta", theta)]:
            z0 = (a[dshift].mean(0) - a[~dshift].mean(0)) / a.std(0)
            dv[f"smd_{name}_abs_max_shifted_labels"] = float(np.abs(z0).max())
            dv[f"smd_{name}_abs_median_shifted_labels"] = float(np.median(np.abs(z0)))
        top = muni.reindex(muni.smd_div_theta.abs().sort_values(ascending=False).index).head(20)
        dv["top20_muni_by_abs_smd_theta"] = top[["CODMUNRES", "births", "deaths", "degree", "smd_div_theta", "smd_div_phi"]].to_dict("records")
        info["divergence_location"] = dv

    # Agreement with what 12_sensitivity_W.py saved for the same variant (0.9 runs reproduce it).
    ref12 = SENS / f"muni_{'queen_refit' if variant == 'queen' else variant}.csv"
    if ref12.exists():
        r = pd.read_csv(ref12, dtype={"CODMUNRES": str})
        fl, fr = (muni["p_gt1.10"] >= 0.95).values, (r["p_gt1.10"] >= 0.95).values
        info["vs_script12_same_variant"] = {
            "max_abs_diff_rr_median": float(np.abs(muni.rr_median.values - r.rr_median.values).max()),
            "max_abs_diff_p110": float(np.abs(muni["p_gt1.10"].values - r["p_gt1.10"].values).max()),
            "flagged_here": int(fl.sum()), "flagged_script12": int(fr.sum()), "flagged_both": int((fl & fr).sum())}
    rq = pd.read_csv(SENS / "muni_queen_refit.csv", dtype={"CODMUNRES": str})
    fl, fq = (muni["p_gt1.10"] >= 0.95).values, (rq["p_gt1.10"] >= 0.95).values
    info["vs_queen_refit_script12"] = {"flagged_here": int(fl.sum()), "flagged_queen": int(fq.sum()), "flagged_both": int((fl & fq).sum()),
                                       "corr_u_mean_vs_log_rr_median_queen": float(np.corrcoef(u_mean, np.log(rq.rr_median.values))[0, 1]),
                                       "sd_log_rr_median_queen": float(np.log(rq.rr_median.values).std())}

    muni.to_csv(DIAG / f"muni_{tag}.csv", index=False)
    np.savez_compressed(DIAG / f"draws_{tag}.npz", diverging=div, mode_ratio=ratio, eigenvalues=lam,
                        **{k: v.astype("float32") for k, v in per_draw.items()})
    (DIAG / f"fit_{tag}.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: v for k, v in info.items() if k not in ("hyper", "divergence_location")}, indent=1), flush=True)
    print(f"{tag}: done {info['elapsed_min']} min, div {info['divergences']}, rhat_hyper {info['rhat_hyper']:.4f}, "
          f"rhat_u {info['rhat_u']:.4f}", flush=True)


def report():
    rows = []
    for f in sorted(DIAG.glob("fit_*.json")):
        i = json.loads(f.read_text())
        r = i["realised"]
        row = {"file": f.stem, "variant": i["variant"], "target_accept": i["target_accept"], "scaling": i["scaling"],
               "divergences": i["divergences"], "div_per_chain": "/".join(map(str, i["divergences_per_chain"])),
               "rhat_hyper": i["rhat_hyper"], "rhat_u": i["rhat_u"],
               "sigma_u": i["hyper"]["sigma_u"]["mean"], "rho": i["hyper"]["rho"]["mean"],
               "ess_bulk_sigma_u": i["hyper"]["sigma_u"]["ess_bulk"], "ess_bulk_rho": i["hyper"]["rho"]["ess_bulk"],
               "sd_u_across_muni": r["sd_u_across_muni"]["mean"], "sd_spatial_across_muni": r["sd_spatial_across_muni"]["mean"],
               "sd_iid_across_muni": r["sd_iid_across_muni"]["mean"], "tau_inv": r["tau_inv_sigma2rho_over_scaling"]["mean"],
               "iid_var": r["iid_var_sigma2_1mrho"]["mean"], "edge_diff_sd_spatial": r["edge_diff_sd_spatial"]["mean"],
               "sd_of_posterior_mean_u": i["sd_of_posterior_mean_u"], "mean_posterior_sd_u": i["mean_posterior_sd_u"],
               "realised_spatial_var": i["realised_spatial_mean_variance_all_modes"],
               "prior_spatial_var_at_post_tau": i["prior_spatial_mean_variance_all_modes_at_posterior_tau"],
               "flagged": i["vs_queen_refit_script12"]["flagged_here"], "flagged_both_queen": i["vs_queen_refit_script12"]["flagged_both"],
               "elapsed_min": i["elapsed_min"]}
        for b in i["spectrum_bands"]:
            row[f"post_over_prior_modes_{b['modes']}"] = b["posterior_over_prior_second_moment_mean"]
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(DIAG / "diagnostics_summary.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 50):
        print(d.round(4).T.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", action="store_true")
    ap.add_argument("--fit", choices=VARIANTS)
    ap.add_argument("--target-accept", type=float, default=0.9)
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--suffix", default="", help="appended to the output tag, e.g. _rerun (same seed, repeat run)")
    a = ap.parse_args()
    if a.graphs:
        graphs()
    elif a.fit:
        fit(a.fit, a.target_accept, a.suffix)
    elif a.report:
        report()
