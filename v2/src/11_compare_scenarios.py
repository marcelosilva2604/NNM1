"""Posterior summaries and comparisons of the missing-input scenarios with the primary fit.

Reads model/idata_bym2_rsr_full.nc (primary) and sensitivity/<scenario>/model/idata_bym2_rsr_full.nc.
Writes outputs/tables/missingness_stress_summary.csv and missingness_stress_coefficients.csv (unrounded values).
Set NNM1_ROOT to the repository's v2 directory and NNM_OUT to redirect the outputs.
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, arviz as az, xarray as xr
from scipy.stats import spearmanr
import nc_io

R = Path(os.environ.get("NNM1_ROOT", Path(__file__).resolve().parents[1]))
OUT = Path(os.environ.get("NNM_OUT", R / "outputs"))


def q(x):
    return dict(zip(["lo95", "median", "hi95"], map(float, np.quantile(x, [0.025, 0.5, 0.975]))))


def summarise(root, name="bym2_rsr_full"):
    d = pd.read_csv(root / "data" / "analysis_muni.csv", dtype={"CODMUNRES": str, "regsaud_id": str})
    post = nc_io.load(root / "model" / f"idata_{name}.nc")["posterior"].to_dataset()
    cols = list(post["beta"].coords["cov"].values)
    beta = post.beta.stack(s=("chain", "draw")).values
    b0 = post.b0.values.ravel()
    u = post.u.stack(s=("chain", "draw")).values
    z = d[cols].to_numpy(float)
    sd = z.std(0)
    z = (z - z.mean(0)) / sd
    rr = np.exp(u)
    expected = d.live_births.to_numpy()[:, None] * np.exp(b0[None, :] + z @ beta)
    departure = expected * (rr - 1)
    p = (rr > 1.1).mean(1)
    flag = p >= 0.95
    mcse = az.mcse(xr.Dataset({"indicator": (post.u.dims, (np.exp(post.u.values) > 1.1).astype(float))}), method="mean")["indicator"].values
    summary = {"n_above": int(flag.sum()), "sd_u": q(u.std(0)), "net": q(departure.sum(0)), "selected": q(departure[flag].sum(0))}
    reference = pd.read_csv(R / "data" / "analysis_muni.csv")
    coef = []
    for j, c in enumerate(cols):
        ref_sd = float(reference[c].std(ddof=0))
        row = {"covariate": c, "scenario_sd": float(sd[j]), "reference_sd": ref_sd, **q(np.exp(beta[j]))}
        row.update({"fixed_sd_" + k: v for k, v in q(np.exp(beta[j] * ref_sd / sd[j])).items()})
        coef.append(row)
    municipal = pd.DataFrame({"CODMUNRES": d.CODMUNRES, "rr_median": np.median(rr, 1), "p": p, "p_mcse": mcse, "flag": flag})
    return summary, municipal, coef


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    primary = summarise(R)
    reference = primary[1].set_index("CODMUNRES")
    rows, coefficients = [], []
    for root in sorted((R / "sensitivity").glob("*")):
        if not (root / "model" / "fit_bym2_rsr_full.json").exists():
            continue
        summary, municipal, coef = summarise(root)
        m = municipal.set_index("CODMUNRES").reindex(reference.index)
        fa, fb = reference.flag, m.flag
        rows.append({"scenario": root.name, "flags": summary["n_above"], "retained": int((fa & fb).sum()),
                     "removed": int((fa & ~fb).sum()), "added": int((~fa & fb).sum()),
                     "spearman_rr": float(spearmanr(reference.rr_median, m.rr_median).statistic),
                     "max_abs_delta_probability": float((m.p - reference.p).abs().max()),
                     **{f"sd_u_{k}": v for k, v in summary["sd_u"].items()},
                     **{f"net_{k}": v for k, v in summary["net"].items()},
                     **{f"selected_{k}": v for k, v in summary["selected"].items()}})
        coefficients += [{"scenario": root.name, "model": "bym2_rsr_full", **c} for c in coef]
    pd.DataFrame(rows).to_csv(OUT / "tables" / "missingness_stress_summary.csv", index=False)
    pd.DataFrame(coefficients).to_csv(OUT / "tables" / "missingness_stress_coefficients.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
