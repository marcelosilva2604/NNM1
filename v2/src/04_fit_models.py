"""Hierarchical spatial models for neonatal mortality beyond structure (artigo 1, v2).

All models share one likelihood: neonatal deaths in municipality i are Poisson with mean
births_i * exp(eta_i), and eta_i = b0 + X_i beta + u_i. The municipality term u_i is the
estimand: the log relative risk beyond what the covariates predict. What changes between
models is the covariate set and the prior on u.

    bym2_rsr_full  X = structure (IDHM income, IDHM education, log population) + maternal
                   composition (mothers < 20 y, >= 35 y, <= 7 y schooling, non-white);
                   u = BYM2 (scaled ICAR spatial + exchangeable) restricted to the orthogonal
                   complement of the covariates (restricted spatial regression). PRIMARY.
    bym2_rsr_struct   same prior, structure only: what composition adds to the map.
    bym2_rsr_full_hr  primary plus an exchangeable health-region effect, for the variance
                   partition (region vs municipality).
    iid_struct     structure only, u exchangeable: the v1 negative binomial fit written as a
                   Poisson-lognormal.
    bym2_*         the same three fits with an unrestricted BYM2 prior, kept as the spatial-
                   confounding sensitivity: the smooth field absorbs the covariates.
    cause_<group>  iid model with the full X, one per action group, outcome = group count.
                   Negative-control contrast: malformations should show the least
                   unexplained heterogeneity if the residual reflects care.

Covariates enter standardised (z-scores over the 5,565 municipalities). Priors are weakly
informative on that scale: beta ~ N(0, 0.5), sigma_u ~ HalfNormal(0.5), rho ~ Beta(1, 1).
The ICAR component is divided by sqrt(scaling factor) so that rho is the share of the
variance of u that is spatially structured (Riebler et al. 2016).

Convention: 95% equal-tailed intervals everywhere.

Outputs (v2/model/):
    idata_<model>.nc        full posterior (b0, beta, sigma_u, rho, u, and gamma_hr where used)
    fit_<model>.json        sampler diagnostics and hyperparameter summaries
    muni_<model>.csv        per-municipality posterior summaries of RR = exp(u)

Run:
    python v2/src/04_fit_models.py [model ...] [--trial]
"""

import json
import sys
import time
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import xarray as xr

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
DATA, OUT = V2 / "data", V2 / "model"
OUT.mkdir(exist_ok=True)

SEED = 20261013
STRUCT = ["idhm_renda", "idhm_educ", "logpop"]
COMPOSITION = ["sh_age_lt20", "sh_age_ge35", "sh_esc_le7", "sh_nonwhite"]
ACTION = ["prenatal", "delivery", "newborn", "malformation"]

ACCESS = ["dist_nicu_100km", "has_nicu", "sh_born_outside", "sh_born_nicu_muni"]
INTERMEDIATE = ["sh_cons_lt7", "sh_preterm", "sh_lbw"]
COMPLETENESS = ["sh_hospital", "sh_miss_esc", "sh_miss_race"]

SPECS = {
    # PRIMARY: restricted spatial regression, structure + composition.
    "bym2_rsr_full": dict(X=STRUCT + COMPOSITION, prior="bym2_rsr", outcome="deaths_total", hr=False),
    "bym2_rsr_struct": dict(X=STRUCT, prior="bym2_rsr", outcome="deaths_total", hr=False),
    "bym2_rsr_full_hr": dict(X=STRUCT + COMPOSITION, prior="bym2_rsr", outcome="deaths_total", hr=True),
    # Unrestricted BYM2: sensitivity for spatial confounding.
    "iid_struct": dict(X=STRUCT, prior="iid", outcome="deaths_total", hr=False),
    "bym2_struct": dict(X=STRUCT, prior="bym2", outcome="deaths_total", hr=False),
    "bym2_full": dict(X=STRUCT + COMPOSITION, prior="bym2", outcome="deaths_total", hr=False),
    "bym2_full_hr": dict(X=STRUCT + COMPOSITION, prior="bym2", outcome="deaths_total", hr=True),
    # Sensitivity fits, exchangeable prior (coefficient questions, not maps):
    #   access         distance to NICU (per 100 km), own NICU, births outside the municipality,
    #                  births in a NICU municipality
    #   intermediate   antenatal visits, preterm, low birth weight: on the pathway, so they
    #                  absorb part of what care does; reported only to show how much
    #   completeness   hospital-birth share and share of births with missing schooling or race/colour:
    #                  registration indicators (descriptive; none measures under-registration)
    "iid_full": dict(X=STRUCT + COMPOSITION, prior="iid", outcome="deaths_total", hr=False),
    "iid_full_access": dict(X=STRUCT + COMPOSITION + ACCESS, prior="iid", outcome="deaths_total", hr=False),
    "iid_full_intermediate": dict(X=STRUCT + COMPOSITION + INTERMEDIATE, prior="iid", outcome="deaths_total", hr=False),
    "iid_full_completeness": dict(X=STRUCT + COMPOSITION + COMPLETENESS, prior="iid", outcome="deaths_total", hr=False),
    **{f"cause_{c}": dict(X=STRUCT + COMPOSITION, prior="iid", outcome=c, hr=False) for c in ACTION},
}
ORDER = ["bym2_rsr_full", "bym2_rsr_struct", "bym2_rsr_full_hr", "iid_struct", "iid_full", "bym2_struct", "bym2_full",
         "bym2_full_hr", "iid_full_access", "iid_full_intermediate", "iid_full_completeness"] + [f"cause_{c}" for c in ACTION]


def load():
    df = pd.read_csv(DATA / "analysis_muni.csv", dtype={"CODMUNRES": str, "regsaud_id": str})
    pairs = pd.read_csv(DATA / "adjacency_pairs.csv")
    scaling = json.loads((DATA / "bym2_scaling.json").read_text())["scaling_factor"]
    n = len(df)
    W = np.zeros((n, n))
    W[pairs.i.values, pairs.j.values] = 1.0
    W[pairs.j.values, pairs.i.values] = 1.0
    return df, W, scaling


def build(df, W, scaling, spec):
    X = df[spec["X"]].values.astype(float)
    Xz = (X - X.mean(0)) / X.std(0)
    n, k = Xz.shape
    y = df[spec["outcome"]].values.astype(int)
    off = np.log(df.live_births.values.astype(float))
    hr_codes, hr_idx = np.unique(df.regsaud_id.values, return_inverse=True)

    with pm.Model(coords={"muni": df.CODMUNRES.values, "cov": spec["X"], "hr": hr_codes}) as model:
        b0 = pm.Normal("b0", np.log(y.sum() / df.live_births.sum()), 1.0)
        beta = pm.Normal("beta", 0.0, 0.5, dims="cov")
        sigma_u = pm.HalfNormal("sigma_u", 0.5)
        theta = pm.Normal("theta", 0.0, 1.0, dims="muni")
        if spec["prior"] in ("bym2", "bym2_rsr"):
            rho = pm.Beta("rho", 1.0, 1.0)
            phi = pm.ICAR("phi", W=W, sigma=1.0, dims="muni")
            conv = pm.math.sqrt(rho / scaling) * phi + pm.math.sqrt(1.0 - rho) * theta
        else:
            conv = theta
        u = sigma_u * conv
        if spec["hr"]:
            sigma_hr = pm.HalfNormal("sigma_hr", 0.5)
            z_hr = pm.Normal("z_hr", 0.0, 1.0, dims="hr")
            gamma = pm.Deterministic("gamma_hr", sigma_hr * z_hr, dims="hr")
            u = u + gamma[hr_idx]
        if spec["prior"] == "bym2_rsr":
            # Restricted spatial regression (Hughes & Haran 2013): the random effect is projected
            # onto the orthogonal complement of the covariate space (intercept included), so the
            # covariates keep the interpretation they have in the non-spatial model and u is, by
            # construction, what structure and composition cannot predict. Without this, the
            # smooth spatial field competes with the smooth covariates and absorbs them
            # (unrestricted fit: income IRR shrinks towards 1, sigma_u doubles, rho -> 0.95).
            Xa = np.column_stack([np.ones(n), Xz])
            M = np.linalg.inv(Xa.T @ Xa)
            u = u - pm.math.dot(Xa, pm.math.dot(M, pm.math.dot(Xa.T, u)))
        pm.Deterministic("u", u, dims="muni")
        eta = b0 + pm.math.dot(Xz, beta) + u
        pm.Poisson("y", mu=pm.math.exp(off + eta), observed=y)
    return model


def summarise(idata, df, name, spec, elapsed, draws, tune, chains):
    post = idata["posterior"]
    u = post["u"].stack(s=("chain", "draw")).values  # n x S
    rr = np.exp(u)
    q = lambda a, p: np.quantile(a, p, axis=1)
    muni = pd.DataFrame({
        "CODMUNRES": df.CODMUNRES.values, "municipio": df.municipio.values, "state": df.state.values,
        "REG": df.REG.values, "deaths": df[spec["outcome"]].values, "births": df.live_births.values,
        "rr_mean": rr.mean(1), "rr_median": q(rr, .5), "rr_lo95": q(rr, .025), "rr_hi95": q(rr, .975),
        "p_rr_gt1": (rr > 1).mean(1), "p_rr_gt1_10": (rr > 1.10).mean(1), "p_rr_gt1_25": (rr > 1.25).mean(1),
        "p_rr_lt0_90": (rr < 0.90).mean(1),
    })
    muni.to_csv(OUT / f"muni_{name}.csv", index=False)

    summ = az.summary(idata, var_names=[v for v in ["b0", "beta", "sigma_u", "rho", "sigma_hr"] if v in post],
                      ci_prob=0.95, round_to=6)
    cols = ["mean", "sd", "eti95_lb", "eti95_ub", "r_hat", "ess_bulk"]
    rhat_u = float(az.rhat(idata, var_names=["u"])["u"].max())
    ess_u = float(az.ess(idata, var_names=["u"])["u"].min())
    info = {
        "model": name, "spec": spec, "n": int(len(df)), "draws": draws, "tune": tune, "chains": chains,
        "elapsed_s": round(elapsed, 1),
        "divergences": int(idata["sample_stats"]["diverging"].sum()),
        "max_rhat_hyper": float(summ["r_hat"].max()), "min_ess_bulk_hyper": float(summ["ess_bulk"].min()),
        "max_rhat_u": rhat_u, "min_ess_bulk_u": ess_u,
        "hyper": summ[cols].to_dict("index"),
    }
    (OUT / f"fit_{name}.json").write_text(json.dumps(info, indent=2, default=str))
    print(json.dumps({k: v for k, v in info.items() if k != "hyper"}, indent=1))
    print(summ[cols].to_string())


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    trial = "--trial" in sys.argv
    names = args or ORDER
    draws, tune, chains = (300, 300, 2) if trial else (1000, 2000, 4)
    df, W, scaling = load()
    print(f"n = {len(df)}, scaling factor = {scaling:.4f}, trial = {trial}", flush=True)
    for name in names:
        spec = SPECS[name]
        path = OUT / f"idata_{name}.nc"
        if path.exists() and not trial:
            print(f"[{name}] cached, skipping", flush=True)
            continue
        t0 = time.time()
        print(f"[{name}] sampling {chains} chains x ({tune} tune + {draws} draws) ...", flush=True)
        with build(df, W, scaling, spec):
            raw = pm.sample(draws=draws, tune=tune, chains=chains, cores=chains, random_seed=SEED,
                            target_accept=0.9, progressbar=False)
        elapsed = time.time() - t0
        print(f"[{name}] done in {elapsed / 60:.1f} min", flush=True)
        # theta and phi are nuisance here; u carries what the paper uses. ArviZ 1.x returns an
        # xarray DataTree, so the subset is rebuilt as a tree before saving.
        ds = raw["posterior"].to_dataset()
        keep = [v for v in ["b0", "beta", "sigma_u", "rho", "sigma_hr", "gamma_hr", "u"] if v in ds]
        ds = ds[keep]
        ds["u"] = ds["u"].astype("float32")
        idata = xr.DataTree.from_dict({"posterior": ds, "sample_stats": raw["sample_stats"].to_dataset()})
        if not trial:
            idata.to_netcdf(path, engine="h5netcdf")
        summarise(idata, df, name if not trial else f"{name}_trial", spec, elapsed, draws, tune, chains)


if __name__ == "__main__":
    main()
