"""Manuscript tables in markdown, built from the result files so no number is typed by hand.

Outputs: v2/outputs/tables/manuscript_tables.md

Run:
    "../artigo 2/.venv/bin/python" v2/src/06_manuscript_tables.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

V2 = Path(__file__).resolve().parents[1]
TAB = V2 / "outputs" / "tables"
N = json.loads((V2 / "outputs" / "NUMBERS.json").read_text())
hyper = pd.read_csv(TAB / "model_hyperparameters_all.csv")
mc = pd.read_csv(TAB / "model_comparison.csv").set_index("model")
df = pd.read_csv(V2 / "data" / "analysis_muni.csv", dtype={"CODMUNRES": str})
m = pd.read_csv(TAB / "municipal_estimates.csv", dtype={"CODMUNRES": str})
unit = pd.read_csv(TAB / "excess_by_unit.csv").set_index("unit")
cause = pd.read_csv(TAB / "cause_negative_control.csv").set_index("cause")

ROWS = [("idhm_renda", "Income sub-index of the municipal HDI"), ("idhm_educ", "Education sub-index of the municipal HDI"),
        ("logpop", "Log resident population"), ("sh_age_lt20", "Mothers younger than 20 years"),
        ("sh_age_ge35", "Mothers aged 35 years or older"), ("sh_esc_le7", "Mothers with 7 or fewer years of schooling"),
        ("sh_nonwhite", "Non-white mothers (black, brown or indigenous)")]
MODELS = [("iid_struct", "Structure, exchangeable"), ("bym2_rsr_struct", "Structure, spatial"), ("bym2_rsr_full", "Structure and composition, spatial (primary)")]


def irr(model, cov):
    r = hyper[(hyper.model == model) & (hyper.covariate == cov)]
    if r.empty:
        return "NA"
    r = r.iloc[0]
    return f"{r.irr:.3f} ({r.lo95:.3f}–{r.hi95:.3f})"


def fmt_ci(a, b, c, d=0):
    return f"{a:,.{d}f} ({b:,.{d}f} to {c:,.{d}f})" if b < 0 or c < 0 else f"{a:,.{d}f} ({b:,.{d}f}–{c:,.{d}f})"


out = []
sd = {k: float(df[k].std()) for k, _ in ROWS}
out.append("**Table 1.** Covariate effects and residual heterogeneity under three specifications. Incidence-rate ratios per standard deviation of the covariate, with 95% credible intervals; the standard deviation is given in parentheses after each covariate. All models are Poisson with a log-births offset; the spatial models carry a BYM2 random effect restricted to the orthogonal complement of the covariates.\n")
out.append("| | " + " | ".join(lab for _, lab in MODELS) + " |")
out.append("|---|" + "---|" * len(MODELS))
for cov, lab in ROWS:
    out.append(f"| {lab} (SD {sd[cov]:.3g}) | " + " | ".join(irr(mod, cov) for mod, _ in MODELS) + " |")
out.append("| SD of log relative risk across municipalities | " + " | ".join(
    f"{mc.loc[mod, 'u_sd_realised_median']:.3f} ({mc.loc[mod, 'u_sd_realised_lo95']:.3f}–{mc.loc[mod, 'u_sd_realised_hi95']:.3f})" for mod, _ in MODELS) + " |")
out.append("| BYM2 spatial share (ρ, prior parameter) | " + " | ".join(
    ("NA" if np.isnan(mc.loc[mod, "rho_median"]) else f"{mc.loc[mod, 'rho_median']:.2f} ({mc.loc[mod, 'rho_lo95']:.2f}–{mc.loc[mod, 'rho_hi95']:.2f})") for mod, _ in MODELS) + " |")
out.append("| Difference in expected log predictive density (WAIC) from the primary model | " + " | ".join(
    f"{mc.loc[mod, 'delta_elpd_vs_primary']:+.0f}" if mod != "bym2_rsr_full" else "0 (reference)" for mod, _ in MODELS) + " |")
out.append("| Municipalities above expectation (P[RR > 1.10] ≥ 0.95) | " + " | ".join(f"{int(mc.loc[mod, 'n_above_10'])}" for mod, _ in MODELS) + " |")
out.append("| Municipalities below expectation (P[RR < 0.90] ≥ 0.95) | " + " | ".join(f"{int(mc.loc[mod, 'n_below_10'])}" for mod, _ in MODELS) + " |")
out.append("")

# Table 2: regions
REG = ["North", "Northeast", "Central-West", "Southeast", "South"]
g = df.groupby("REG").agg(n=("CODMUNRES", "size"), deaths=("deaths_total", "sum"), births=("live_births", "sum"))
out.append("**Table 2.** Neonatal mortality and the residual relative risk by macro-region, 2014–2024. Excess deaths are the posterior median of Σ E_i(RR_i − 1) over the municipalities in each set, with 95% credible intervals; the net excess sums over all municipalities of the region and can be negative.\n")
out.append("| Region | Municipalities | Neonatal deaths | Rate per 1,000 live births | Above expectation, n | Below expectation, n | Excess deaths in municipalities above expectation | Net excess deaths, all municipalities |")
out.append("|---|---|---|---|---|---|---|---|")
for r in REG + ["Brazil"]:
    if r == "Brazil":
        n_, d_, b_ = len(df), int(df.deaths_total.sum()), int(df.live_births.sum())
    else:
        n_, d_, b_ = int(g.loc[r, "n"]), int(g.loc[r, "deaths"]), int(g.loc[r, "births"])
    u = unit.loc[r]
    out.append(f"| {r} | {n_:,} | {d_:,} | {d_ / b_ * 1000:.2f} | {int(u.n_above)} | {int(u.n_below)} | "
               f"{fmt_ci(u.excess_in_above_median, u.excess_in_above_lo95, u.excess_in_above_hi95)} | "
               f"{fmt_ci(u.net_median, u.net_lo95, u.net_hi95)} |")
out.append("")

# Table 3: top 12 by excess deaths among the above set
t = pd.read_csv(TAB / "top20_excess_deaths.csv").head(12)
out.append("**Table 3.** The twelve municipalities with the largest excess deaths among those above expectation with at least 95% probability. Expected deaths come from structure and maternal composition; the relative risk is the residual beyond them. Hospital births is the share of the municipality's live births that took place in a hospital, a record-quality and access indicator.\n")
out.append("| Municipality (state) | Neonatal deaths | Expected deaths | Relative risk (95% CrI) | Excess deaths (95% CrI) | Hospital births (%) |")
out.append("|---|---|---|---|---|---|")
for _, r in t.iterrows():
    out.append(f"| {r.municipio} ({r.state}) | {int(r.deaths):,} | {r.expected_median:,.0f} | {r.rr_median:.2f} ({r.rr_lo95:.2f}–{r.rr_hi95:.2f}) | "
               f"{r.excess_median:,.0f} ({r.excess_lo95:,.0f}–{r.excess_hi95:,.0f}) | {r.sh_hospital * 100:.0f} |")
out.append("")

# Supplementary pieces
out.append("**Table S1.** All fitted models: residual heterogeneity, spatial share, health-region standard deviation, fit and classification counts.\n")
out.append("| Model | Covariates | Prior on u | SD of log RR (95% CrI) | ρ | Health-region effect, SD on the log scale | Δelpd (WAIC) vs primary | Above | Below | Spearman with primary |")
out.append("|---|---|---|---|---|---|---|---|---|---|")
DESC = {"iid_struct": ("structure", "exchangeable"), "iid_full": ("structure + composition", "exchangeable"),
        "bym2_rsr_struct": ("structure", "BYM2, restricted"), "bym2_rsr_full": ("structure + composition", "BYM2, restricted (primary)"),
        "bym2_rsr_full_hr": ("structure + composition", "BYM2, restricted, + health region"),
        "bym2_struct": ("structure", "BYM2, unrestricted"), "bym2_full": ("structure + composition", "BYM2, unrestricted"),
        "bym2_full_hr": ("structure + composition", "BYM2, unrestricted, + health region"),
        "iid_full_access": ("+ access to care", "exchangeable"), "iid_full_intermediate": ("+ antenatal visits, preterm, low birth weight", "exchangeable"),
        "iid_full_completeness": ("+ registration indicators (hospital births, missing schooling and race/colour)", "exchangeable")}
for mod, (cv, pr) in DESC.items():
    r = mc.loc[mod]
    rho = "NA" if np.isnan(r.rho_median) else f"{r.rho_median:.2f} ({r.rho_lo95:.2f}–{r.rho_hi95:.2f})"
    icc = "NA" if np.isnan(r.get("sigma_hr_median", np.nan)) else f"{r.sigma_hr_median:.3f} ({r.sigma_hr_lo95:.3f}–{r.sigma_hr_hi95:.3f})"
    out.append(f"| {mod} | {cv} | {pr} | {r.u_sd_realised_median:.3f} ({r.u_sd_realised_lo95:.3f}–{r.u_sd_realised_hi95:.3f}) | {rho} | {icc} | "
               f"{r.delta_elpd_vs_primary:+.0f} | {int(r.n_above_10)} | {int(r.n_below_10)} | {r.spearman_u_vs_primary:.2f} |")
out.append("")
out.append("**Table S2.** Sensitivity covariates (exchangeable models, incidence-rate ratio per SD with 95% CrI).\n")
out.append("| Model | Covariate (SD) | IRR (95% CrI) |")
out.append("|---|---|---|")
for mod in ["iid_full_access", "iid_full_intermediate", "iid_full_completeness"]:
    s = hyper[(hyper.model == mod) & (hyper.covariate != "") & ~hyper.term.isin(["sigma_u", "rho", "sigma_hr", "diagnostics"])]
    for _, r in s.iterrows():
        out.append(f"| {mod} | {r.term} ({r.unit.replace('per SD (', '').rstrip(')')}) | {r.irr:.3f} ({r.lo95:.3f}–{r.hi95:.3f}) |")
out.append("")
out.append("**Table S3.** Residual heterogeneity by cause group (exchangeable model with structure and composition, one fit per group) and the cause-specific relative risk in the municipalities classified above expectation on all causes.\n")
out.append("| Cause group | Deaths | SD of log RR (95% CrI) | Median cause-specific RR in the 'above' municipalities | Median in the others | Municipalities above expectation on this cause |")
out.append("|---|---|---|---|---|---|")
LAB = {"prenatal": "Antenatal-care causes", "delivery": "Delivery-care causes", "newborn": "Newborn-care causes", "malformation": "Congenital malformations"}
for c in ["prenatal", "delivery", "newborn", "malformation"]:
    r = cause.loc[c]
    out.append(f"| {LAB[c]} | {int(r.deaths):,} | {r.sigma_u_median:.3f} ({r.sigma_u_lo95:.3f}–{r.sigma_u_hi95:.3f}) | {r.rr_median_in_above_set:.2f} | {r.rr_median_in_others:.2f} | {int(r.n_above_10)} |")
out.append("")
det = pd.read_csv(TAB / "detectability.csv")
out.append("**Table S4.** Minimum residual relative risk detectable with 80% power (one-sided α = 0.05) from sampling variation alone, by expected deaths, and the number of municipalities at or above each expected count.\n")
out.append("| Expected deaths | Minimum detectable RR | Municipalities with at least this expectation | Share of all neonatal deaths in them |")
out.append("|---|---|---|---|")
for _, r in det.iterrows():
    out.append(f"| {int(r.expected_deaths):,} | {r.mde_rr_poisson:.2f} | {int(r.municipalities_at_or_above):,} | {r.deaths_share_at_or_above * 100:.0f}% |")
out.append("")
p = N["ppc"]
out.append(f"**Posterior predictive check (primary model).** {p['share_outside_95pp'] * 100:.1f}% of municipalities fall outside their 95% posterior predictive interval; "
           f"share of municipalities with zero deaths {p['share_zero_obs'] * 100:.1f}% observed versus {p['share_zero_pp'] * 100:.1f}% predicted; "
           f"largest count {p['max_obs']:,} observed versus {p['max_pp_median']:,.0f} predicted (median); variance-to-mean ratio {p['var_mean_ratio_obs']:,.0f} versus {p['var_mean_ratio_pp_median']:,.0f}.")
(TAB / "manuscript_tables.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
print("\nregional NMR check:", {r: round(g.loc[r, 'deaths'] / g.loc[r, 'births'] * 1000, 2) for r in REG})
print("composition SDs:", {k: round(v, 4) for k, v in sd.items()})
print("v1 ids in analysis:", int(df.v1_included.sum()))
