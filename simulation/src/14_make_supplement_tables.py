"""Build the supplementary tables (markdown) from the grid outputs. Read-only on inputs; writes outputs/supp_tables.md."""
import glob
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
O = HERE.parent / "outputs"
d = pd.read_csv(O / "grid_summary.csv")
m = pd.read_csv(O / "grid_replicate_meta.csv")
cal = pd.read_csv(O / "grid_calibration.csv")
KS = [0.5, 1.0, 1.5, 2.0]
pct = lambda x: "NA" if pd.isna(x) else f"{100 * x:.1f}%"
out = []


def md(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.values]
    return "\n".join(lines)


# eTable 1: calibration
t = cal[cal.k > 0].copy()
t["bin"] = t.apply(lambda r: f"{r.p_lo:.2f}-{r.p_hi:.2f}", axis=1)
a = t.pivot(index="bin", columns="k", values="share_truly_gt110").map(pct)
b = t.pivot(index="bin", columns="k", values="n")
rows = []
for bin_ in a.index:
    rows.append([bin_] + [f"{a.loc[bin_, k]} (n = {int(b.loc[bin_, k]):,})" for k in KS])
out.append("**eTable 1.** Share of municipality-replicates with true RR > 1.10 by bin of estimated P(RR > 1.10), with the number of municipality-replicates in each bin.\n\n" + md(pd.DataFrame(rows, columns=["Estimated probability"] + [f"k = {k}" for k in KS])))

# eTable 2: rule family and comparators, all k
rules = [("post_tau1.05_c0.95", "P(RR > 1.05) >= 0.95"), ("post_tau1.10_c0.99", "P(RR > 1.10) >= 0.99"), ("post_tau1.10_c0.95", "P(RR > 1.10) >= 0.95 (primary)"),
         ("post_tau1.10_c0.90", "P(RR > 1.10) >= 0.90"), ("post_tau1.10_c0.80", "P(RR > 1.10) >= 0.80"), ("post_tau1.20_c0.95", "P(RR > 1.20) >= 0.95"),
         ("C1_median_gt1.10", "C1: posterior median RR > 1.10"), ("C2_poisson_p05", "C2: Poisson tail test, p < 0.05"), ("C3_poisson_BH05", "C3: Poisson tail test, BH 5%")]
rows = []
for r, lab in rules:
    x = d[(d.rule == r) & (d.stratum == "ALL")].set_index("k")
    rows.append([lab] + [f"{pct(x.loc[k, 'sens_pooled'])} / {pct(x.loc[k, 'false_flag_pooled'])} / {int(x.loc[k, 'n_flag_total']):,}" for k in KS])
out.append("**eTable 2.** Rules and comparators by k, shown as sensitivity / false-flag proportion / flagged municipality-replicates (pooled; target true RR above the rule's own tau, which is 1.10 except in the first and sixth rows). Under k = 0, C2 flagged 3,159 municipality-replicates in 40 replicates; the other rules flagged none.\n\n" + md(pd.DataFrame(rows, columns=["Rule"] + [f"k = {k}" for k in KS])))

# eTable 3: coverage by stratum
rows = []
reps = [json.loads(Path(f).read_text()) for f in glob.glob(str(O / "grid" / "rep_k*.json"))]
agg = {}
for r in reps:
    for lab, v in r["coverage_strata"].items():
        z = agg.setdefault((r["k"], lab), [0, 0]); z[0] += v["n"]; z[1] += v["n_cov"]
for lab in ["<20", "20-49", "50-199", ">=200"]:
    rows.append([lab] + [pct(agg[(k, lab)][1] / agg[(k, lab)][0]) for k in [0.0] + KS])
# coverage among flagged municipality-replicates (primary rule): pooled and mean of replicate-level values
fl = {}
for r in reps:
    nf = r["rules"]["post_tau1.10_c0.95"]["n_flag"]
    if nf and r["coverage_flagged"] is not None:
        z = fl.setdefault(r["k"], [0, 0.0, []]); z[0] += nf; z[1] += r["coverage_flagged"] * nf; z[2].append(r["coverage_flagged"])
rows.append(["Among flagged, pooled", "NA"] + [pct(fl[k][1] / fl[k][0]) for k in KS])
rows.append(["Among flagged, mean of replicate-level values", "NA"] + [pct(sum(fl[k][2]) / len(fl[k][2])) for k in KS])
out.append("**eTable 3.** Pooled coverage of the 95% credible interval for the residual RR, by generating expected deaths, and among flagged municipality-replicates under the primary rule.\n\n" + md(pd.DataFrame(rows, columns=["Stratum or subset", "k = 0"] + [f"k = {k}" for k in KS])))

# eTable 4: convergence by k
rows = []
for k in [0.0] + KS:
    x = m[m.k == k]
    hy = x.rhat_hyper.notna()
    rows.append([k, len(x), int((x.divergences > 0).sum()), int(x.divergences.sum()), f"{int((x.rhat_hyper > 1.10).sum())} of {int(hy.sum())}",
                 f"{x.rhat_hyper.max():.3f}" if hy.any() else "NA", f"{int((x.rhat_u > 1.10).sum())} of {int(x.rhat_u.notna().sum())}", f"{x.rhat_u.max():.3f}",
                 f"{int((x.rhat_max > 1.10).sum())} of {len(x)}", f"{x.rhat_max.max():.3f}"])
out.append("**eTable 4.** Convergence diagnostics by scenario. The hyperparameter and municipal-effect R-hat were recorded only for replicates run after the plan amendment.\n\n" + md(pd.DataFrame(rows, columns=["k", "Replicates", "With divergences", "Divergences", "Hyperparameter R-hat > 1.10", "Max hyperparameter R-hat", "Municipal-effect R-hat > 1.10", "Max municipal-effect R-hat", "Overall R-hat > 1.10", "Max overall R-hat"])))

# eTable 5: generating heterogeneity and target prevalence
rows = []
for k in KS:
    x = m[m.k == k]; y = d[(d.rule == "post_tau1.10_c0.95") & (d.stratum == "ALL") & (d.k == k)].iloc[0]
    rows.append([k, f"{x.sd_u_true.mean():.3f} ({x.sd_u_true.min():.3f}-{x.sd_u_true.max():.3f})", pct(y.n_true_total / (100 * 5565))]
                + [pct(d[(d.rule == "post_tau1.10_c0.95") & (d.stratum == s) & (d.k == k)].iloc[0].n_true_total / (100 * n)) for s, n in [("<20", 0), ("20-49", 0), ("50-199", 0), (">=200", 0)]][:0])
out.append("**eTable 5.** Standard deviation of the generating log RR across municipalities (mean and range over replicates) and share of municipality-replicates with true RR > 1.10.\n\n" + md(pd.DataFrame(rows, columns=["k", "SD of generating log RR", "True RR > 1.10"])))
(O / "supp_tables.md").write_text("\n\n".join(out) + "\n")
print("\n\n".join(out))
