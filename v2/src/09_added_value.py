"""What the model adds to routine outputs: comparison of who is flagged.

Compares the municipalities classified above expectation by the primary model (P[RR > 1.10] >= 0.95)
with (a) two structure-only models, non-spatial with an exchangeable random effect (what a standard
expected-count smoothing gives) and spatial with the same restricted BYM2 prior, and (b) the municipalities
with the highest crude neonatal mortality rate among those with at least 50 deaths. Also reports the model-based expected number of flagged
municipalities whose true relative risk is at or below 1.10, sum(1 - P[RR > 1.10]) over the flagged set.

Needs model/muni_iid_struct.csv and model/muni_bym2_rsr_struct.csv (written by 04) and outputs/tables/municipal_estimates.csv (05).
Output: outputs/added_value.json

Run:
    python src/09_added_value.py
"""

import json
import os as _os
from pathlib import Path

import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = Path(_os.environ.get("NNM1_ROOT", HERE.parent))  # repository root; override with NNM1_ROOT

P = pd.read_csv(ROOT / "outputs" / "tables" / "municipal_estimates.csv", dtype={"CODMUNRES": str}).set_index("CODMUNRES")
SP = pd.read_csv(ROOT / "model" / "muni_bym2_rsr_struct.csv", dtype={"CODMUNRES": str}).set_index("CODMUNRES")
C = pd.read_csv(ROOT / "model" / "muni_iid_struct.csv", dtype={"CODMUNRES": str}).set_index("CODMUNRES")
fp = set(P.index[P["class"] == "above"])
fc = set(C.index[C.p_rr_gt1_10 >= 0.95])
fs = set(SP.index[SP.p_rr_gt1_10 >= 0.95])
big = P[P.deaths >= 50].copy()
big["nmr"] = big.deaths / big.births * 1000
k = len(fp)
top = set(big.sort_values("nmr", ascending=False).index[:k])
above = P.loc[list(fp)]
out = {
    "primary_flagged": len(fp), "conventional_flagged": len(fc), "in_common": len(fp & fc),
    "primary_only": len(fp - fc), "conventional_only": len(fc - fp),
    "primary_only_by_region": P.loc[list(fp - fc), "REG"].value_counts().to_dict(),
    "spearman_rr_primary_vs_conventional": float(P.rr_median.corr(C.rr_median.reindex(P.index), method="spearman")),
    "spatial_struct_flagged": len(fs), "spatial_struct_in_common": len(fp & fs), "primary_only_vs_spatial_struct": len(fp - fs),
    "spatial_struct_only": len(fs - fp),
    "crude_nmr_top_k_k": k, "crude_nmr_top_k_in_common_with_primary": len(top & fp),
    "flagged_with_ge50_deaths": int(len(fp & set(big.index))),
    "flagged_ge50_in_top_k_crude": len(top & fp),
    "top67_crude_ge50_overlap_with_flagged_ge50": len(set(big.sort_values("nmr", ascending=False).index[:len(fp & set(big.index))]) & fp),
    "expected_false_flags_model_based": float((1 - above.p_rr_gt1_10).sum()),
    "expected_false_flag_share_model_based": float((1 - above.p_rr_gt1_10).mean()),
}
(ROOT / "outputs" / "added_value.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=1))
