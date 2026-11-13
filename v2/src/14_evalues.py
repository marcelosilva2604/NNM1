"""E-values (VanderWeele & Ding 2017) for the primary model (bym2_rsr_full), v2 outputs.
Reads v2/outputs/tables/model_hyperparameters_all.csv (unrounded posterior IRR and 95% CrI).
E = RR + sqrt(RR*(RR-1)) for RR >= 1; ratios below 1 are inverted first.
For the limit of the interval closer to 1: E = 1 if the interval includes 1."""
import math, pandas as pd, sys
from pathlib import Path
def ev(rr):
    rr = 1/rr if rr < 1 else rr
    return rr + math.sqrt(rr*(rr-1))
def ev_ci(rr, lo, hi):
    if lo <= 1 <= hi: return 1.0
    lim = lo if rr > 1 else hi   # limit closer to the null
    return ev(lim)
# sanity checks of the implementation
assert abs(ev(1.0)-1.0)<1e-12 and abs(ev(2.0)-(2+math.sqrt(2)))<1e-12 and abs(ev(0.5)-ev(2.0))<1e-12
d = pd.read_csv(Path(__file__).resolve().parents[1] / 'outputs/tables/model_hyperparameters_all.csv')
p = d[(d.model=='bym2_rsr_full') & d.irr.notna() & d.covariate.isin(['idhm_renda','idhm_educ','logpop','sh_age_lt20','sh_age_ge35','sh_esc_le7','sh_nonwhite'])]
rows=[]
for _,r in p.iterrows():
    rows.append((r.term, r.irr, r.lo95, r.hi95, ev(r.irr), ev_ci(r.irr,r.lo95,r.hi95)))
out = pd.DataFrame(rows, columns=['term','irr','lo95','hi95','E_point','E_ci_limit'])
pd.set_option('display.width',200)
print(out.round(4).to_string(index=False))
print('\nResidual RR thresholds:')
for t in (1.10,1.25): print(t, round(ev(t),3))
out.to_csv(Path(__file__).resolve().parents[1] / 'outputs/tables/evalues_primary_model.csv', index=False)
