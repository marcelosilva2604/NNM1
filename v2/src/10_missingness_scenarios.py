"""Missing-input scenarios for maternal schooling and race/skin colour (supplement, eMethods 3).

For each category, N = registered births, K = births with a known value, C = known births in the category, M = N - K.
Extreme scenarios use (C + qM)/N with q = 0 or 1 in all four combinations. Intermediate scenarios use
q = w*p/(1 - p + w*p) with p = C/K and w = 0.5 or 2. These are deterministic input scenarios, not imputations and not
bounds on any model output. Each scenario writes its own analysis table under sensitivity/<name>/ and refits the
primary specification with 04_fit_models.py (set NNM_DRAWS and NNM_TUNE for the sampling settings).

Run:  python src/10_missingness_scenarios.py [--intermediate] [--no-fit]
"""
import json, os, shutil, subprocess, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
V2 = Path(os.environ.get("NNM1_ROOT", HERE.parent))
EXTREME = [(0, 0), (0, 1), (1, 0), (1, 1)]
ODDS = [(0.5, 0.5), (0.5, 2), (2, 0.5), (2, 2)]


def scenario_table(base, schooling, race, odds):
    d = base.copy()
    for out, num, den, mult in [("sh_esc_le7", "esc_le7", "esc_known", schooling),
                                ("sh_nonwhite", "race_nonwhite", "race_known", race)]:
        p = d[num] / d[den]
        q = mult * p / (1 - p + mult * p) if odds else mult
        d[out] = (d[num] + q * (d.births - d[den])) / d.births
    return d


def main(intermediate=False, fit=True):
    base = pd.read_csv(V2 / "data" / "analysis_muni.csv", dtype={"CODMUNRES": str, "regsaud_id": str})
    for schooling, race in (ODDS if intermediate else EXTREME):
        tag = lambda x: str(x).replace(".", "p")
        name = f"odds_e{tag(schooling)}_r{tag(race)}" if intermediate else f"extreme_e{schooling}_r{race}"
        root = V2 / "sensitivity" / name
        (root / "data").mkdir(parents=True, exist_ok=True)
        (root / "model").mkdir(exist_ok=True)
        scenario_table(base, schooling, race, intermediate).to_csv(root / "data" / "analysis_muni.csv", index=False)
        for f in ("adjacency_pairs.csv", "bym2_scaling.json"):
            shutil.copy2(V2 / "data" / f, root / "data" / f)
        (root / "scenario.json").write_text(json.dumps({
            "type": "odds_multiplier" if intermediate else "extreme", "schooling": schooling, "race": race,
            "description": "Among records with a missing value, the category probability follows the stated extreme or odds "
                           "multiplier relative to known records within each municipality. Both proportions use all registered "
                           "births as denominator. Deterministic input sensitivity, not multiple imputation."}, indent=2))
        if fit:
            env = dict(os.environ, NNM1_ROOT=str(root), NNM_DRAWS=os.environ.get("NNM_DRAWS", "4000"),
                       NNM_TUNE=os.environ.get("NNM_TUNE", "4000"))
            subprocess.run([sys.executable, str(HERE / "04_fit_models.py"), "bym2_rsr_full"], env=env, check=True)


if __name__ == "__main__":
    main("--intermediate" in sys.argv, "--no-fit" not in sys.argv)
