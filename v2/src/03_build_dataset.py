"""Assemble the v2 municipality table.

One row per municipality of residence, pooled 2014-2024, with the same inclusion rule
as the earlier non-spatial analysis (live births > 0, IDHM and population available), which
gave 5,477 municipalities; the analysis table has 5,565 (88 were not in that analysis), plus the variables v1 did not have:

    maternal composition       shares of mothers < 20 y, >= 35 y, <= 7 y schooling, non-white
    place of birth             shares born in hospital, outside the municipality, in a
                               municipality holding NICU beds that year
    intermediates              < 7 antenatal visits, preterm, low birth weight (descriptive
                               and sensitivity only: they sit on the pathway from care to death)
    registration indicators    share of hospital births, share of births with missing maternal
                               schooling or race/colour (descriptive indicators of record quality;
                               none measures under-registration of births or deaths)
    geography                  state, IBGE immediate region, SUS health region

Access proxies from v1 are recomputed unchanged (distance to the nearest NICU
municipality, NICU beds and primary-care units per 1,000 births), and reported on a
readable scale downstream (per 100 km).

Outputs (v2/data/):
    analysis_muni.csv
    dataset_report.txt

Run:
    python v2/src/03_build_dataset.py
"""

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
A1 = V2.parent
A2REF = Path(__import__("os").environ.get("NNM1_REF", A1 / "data" / "ref"))  # health-region tables; not distributed
PROC = A1 / "data" / "processed"
RAW = A1 / "data" / "raw"

STATE_NAMES = {"11": "RO", "12": "AC", "13": "AM", "14": "RR", "15": "PA", "16": "AP", "17": "TO",
               "21": "MA", "22": "PI", "23": "CE", "24": "RN", "25": "PB", "26": "PE", "27": "AL",
               "28": "SE", "29": "BA", "31": "MG", "32": "ES", "33": "RJ", "35": "SP", "41": "PR",
               "42": "SC", "43": "RS", "50": "MS", "51": "MT", "52": "GO", "53": "DF"}
REGION = {"1": "North", "2": "Northeast", "3": "Southeast", "4": "South", "5": "Central-West"}


def nicu_distance(cn):
    """Haversine distance (km) from each municipality to the nearest one holding NICU beds in the last CNES year."""
    la, lo = defaultdict(list), defaultdict(list)
    with open(RAW / "cnes_estabelecimentos.csv", encoding="latin-1") as fh:
        r = csv.reader(fh, delimiter=";")
        h = [x.strip().strip('"').upper() for x in next(r)]
        ix = {x: i for i, x in enumerate(h)}
        for row in r:
            if len(row) <= ix["NU_LONGITUDE"]:
                continue
            try:
                a = float(row[ix["NU_LATITUDE"]]); o = float(row[ix["NU_LONGITUDE"]])
            except ValueError:
                continue
            if -34 <= a <= 6 and -74 <= o <= -34:
                ib = row[ix["CO_IBGE"]].strip().strip('"')[:6]
                la[ib].append(a); lo[ib].append(o)
    cent = {ib: (np.median(la[ib]), np.median(lo[ib])) for ib in la}
    last = cn.sort_values("ANO").groupby("CODMUNRES").last()
    nc = np.array([cent[m] for m in set(last.index[last["uti_neo_leitos"] > 0]) if m in cent])

    def hav(a1, o1, A, O):
        R = 6371; p = np.pi / 180; d1 = (A - a1) * p; d2 = (O - o1) * p
        return 2 * R * np.arcsin(np.sqrt(np.sin(d1 / 2) ** 2 + np.cos(a1 * p) * np.cos(A * p) * np.sin(d2 / 2) ** 2))

    def dist(m):
        if m not in cent:
            return np.nan
        a, o = cent[m]
        return float(hav(a, o, nc[:, 0], nc[:, 1]).min())

    return dist


def main():
    # v1 built the municipality table from the deaths file, so the 112 municipalities with
    # live births and no neonatal death in the decade were left out (selection on the
    # outcome). Here the universe is every municipality with a live birth; deaths are zero
    # where SIM has none.
    deaths = pd.read_csv(PROC / "muni_base_2014_2024.csv", dtype={"CODMUNRES": str}).drop(columns=["live_births", "nmr"])
    births_all = pd.read_csv(PROC / "nascidos_muni_ano.csv", dtype={"CODMUNRES": str})
    births_all["CODMUNRES"] = births_all.CODMUNRES.str[:6]
    births_all = births_all.groupby("CODMUNRES", as_index=False).nascimentos.sum().rename(columns={"nascimentos": "live_births"})
    base = births_all.merge(deaths, on="CODMUNRES", how="outer")
    for c in ["prenatal", "delivery", "newborn", "malformation", "illdef", "deaths_total"]:
        base[c] = base[c].fillna(0).astype(int)
    base["live_births"] = base.live_births.fillna(0).astype(int)
    n_zero_death = int(((base.deaths_total == 0) & (base.live_births > 0)).sum())
    base["nmr"] = np.where(base.live_births > 0, base.deaths_total / base.live_births.replace(0, np.nan) * 1000, np.nan)
    idhm = pd.read_csv(PROC / "idhm_muni.csv", dtype={"CODMUNRES": str})
    pop = (pd.read_csv(PROC / "populacao_muni_ano.csv", dtype={"CODMUNRES": str})
             .groupby("CODMUNRES").populacao.mean().rename("pop_mean"))
    cn = pd.read_csv(PROC / "cnes_muni_ano.csv", dtype={"CODMUNRES": str})
    cn["CODMUNRES"] = cn.CODMUNRES.str[:6]
    acc = cn.groupby("CODMUNRES").agg(ubs_mean=("ubs", "mean"), nicu_beds_mean=("uti_neo_leitos", "mean"),
                                      nicu_any_year=("uti_neo_leitos", lambda s: int((s > 0).any())))
    dist = nicu_distance(cn)

    mat = pd.read_csv(V2 / "data" / "maternal_composition_muni.csv", dtype={"CODMUNRES": str})
    rs = pd.read_csv(A2REF / "muni_regiao_saude.csv", dtype={"cod6": str, "regsaud_id": str})
    ri = pd.read_csv(A2REF / "muni_regiao_imediata.csv", dtype={"CODMUNRES": str, "rgi_id": str})
    v1 = pd.read_csv(A1 / "outputs" / "tables" / "excess_mortality.csv", dtype={"CODMUNRES": str})
    v1 = v1[["CODMUNRES", "expA", "SMR_A", "SMR_A_eb", "sig_excess", "sig_deficit"]].rename(
        columns={"expA": "v1_expected", "SMR_A": "v1_smr", "SMR_A_eb": "v1_smr_eb",
                 "sig_excess": "v1_sig_excess", "sig_deficit": "v1_sig_deficit"})

    df = (base.merge(idhm[["CODMUNRES", "municipio", "uf", "idhm", "idhm_renda", "idhm_long", "idhm_educ"]],
                     on="CODMUNRES", how="left")
              .merge(pop, on="CODMUNRES", how="left")
              .merge(acc, on="CODMUNRES", how="left")
              .merge(mat, on="CODMUNRES", how="left")
              .merge(rs[["cod6", "regsaud_id", "regsaud_nome"]], left_on="CODMUNRES", right_on="cod6", how="left")
              .merge(ri[["CODMUNRES", "rgi_id", "rgi_nome"]], on="CODMUNRES", how="left")
              .merge(v1, on="CODMUNRES", how="left"))
    df = df.drop(columns=["cod6"])
    df["REG"] = df.CODMUNRES.str[0].map(REGION)
    df["state"] = df.CODMUNRES.str[:2].map(STATE_NAMES)
    df["dist_nicu_km"] = df.CODMUNRES.map(dist)

    included = (df.live_births > 0) & df.idhm_renda.notna() & df.pop_mean.notna()
    excluded = df[~included]
    df = df[included].copy().sort_values("CODMUNRES").reset_index(drop=True)
    df["v1_sig_excess"] = df.v1_sig_excess.fillna(False).astype(bool)
    df["v1_sig_deficit"] = df.v1_sig_deficit.fillna(False).astype(bool)
    df["v1_included"] = df.v1_expected.notna()

    df["logbirths"] = np.log(df.live_births)
    df["logpop"] = np.log(df.pop_mean)
    df["dist_nicu_km"] = df.dist_nicu_km.fillna(df.dist_nicu_km.median())
    df["dist_nicu_100km"] = df.dist_nicu_km / 100.0
    df["nicu_per1k"] = df.nicu_beds_mean.fillna(0) / df.live_births * 1000
    df["ubs_per1k"] = df.ubs_mean.fillna(0) / df.live_births * 1000
    df["has_nicu"] = (df.nicu_beds_mean.fillna(0) > 0).astype(int)

    def share(num, den):
        return np.where(df[den] > 0, df[num] / df[den].replace(0, np.nan), np.nan)

    df["sh_age_lt20"] = share("age_lt20", "age_known")
    df["sh_age_ge35"] = share("age_ge35", "age_known")
    df["sh_esc_le7"] = share("esc_le7", "esc_known")
    df["sh_esc_le3"] = share("esc_le3", "esc_known")
    df["sh_nonwhite"] = share("race_nonwhite", "race_known")
    df["sh_hospital"] = share("loc_hospital", "loc_known")
    df["sh_born_outside"] = share("born_outside", "nasc_known")
    df["sh_born_nicu_muni"] = share("born_in_nicu_muni", "nasc_known")
    df["sh_cons_lt7"] = share("cons_lt7", "cons_known")
    df["sh_preterm"] = share("preterm", "gest_known")
    df["sh_lbw"] = share("lbw", "peso_known")
    df["sh_miss_esc"] = 1 - share("esc_known", "births")
    df["sh_miss_race"] = 1 - share("race_known", "births")
    df["births_diff_sinasc"] = df.births - df.live_births  # SINASC files vs v1 denominator
    df["illdef_share"] = df.illdef / (df.deaths_total + df.illdef).replace(0, np.nan)

    # Shares without a known denominator are set to the national median so that no covariate carries a missing value.
    for c in ["sh_age_lt20", "sh_age_ge35", "sh_esc_le7", "sh_nonwhite", "sh_hospital", "sh_born_outside",
              "sh_born_nicu_muni", "sh_cons_lt7", "sh_preterm", "sh_lbw"]:
        miss = int(df[c].isna().sum())
        if miss:
            df[c] = df[c].fillna(df[c].median())
            print(f"{c}: {miss} municipalities without a defined value, set to the median")

    df.to_csv(V2 / "data" / "analysis_muni.csv", index=False)

    births_national = pd.read_csv(PROC / "nascidos_muni_ano.csv", dtype={"CODMUNRES": str}).nascimentos.sum()
    deaths_national = base.deaths_total.sum()
    rel = df.deaths_total >= 50
    cov = ["idhm_renda", "idhm_educ", "logpop", "sh_age_lt20", "sh_age_ge35", "sh_esc_le7", "sh_nonwhite",
           "sh_hospital", "sh_born_outside", "sh_born_nicu_muni", "sh_cons_lt7", "sh_preterm", "sh_lbw",
           "dist_nicu_100km", "sh_miss_esc", "sh_miss_race"]
    corr = df.loc[rel, cov + ["v1_smr"]].corr(method="spearman")["v1_smr"].drop("v1_smr")
    lines = [
        "v2 analysis table", "=" * 60,
        f"national: {deaths_national:,} neonatal deaths / {births_national:,} live births (SINASC complete) "
        f"= {deaths_national / births_national * 1000:.2f} per 1,000",
        f"municipalities with births: {int((base.live_births > 0).sum()):,}; with births and no neonatal death: {n_zero_death}",
        f"included municipalities : {len(df):,}  ({int(df.deaths_total.sum()):,} deaths, {int(df.live_births.sum()):,} births; "
        f"{int((df.deaths_total == 0).sum())} with zero deaths; {int((~df.v1_included).sum())} not in v1)",
        f"excluded                : {len(excluded)} ({excluded.CODMUNRES.tolist()}; births {int(excluded.live_births.sum()):,}, deaths {int(excluded.deaths_total.sum()):,})",
        f"with >= 50 deaths       : {int(rel.sum())}",
        f"births, linked files vs v1 denominator: total diff {int(df.births_diff_sinasc.sum()):,}; "
        f"max |diff| {int(df.births_diff_sinasc.abs().max()):,}",
        f"registration indicators: median missing schooling {df.sh_miss_esc.median():.3f}, missing race/colour {df.sh_miss_race.median():.3f}, "
        f"hospital births {df.sh_hospital.median():.3f}",
        "",
        "covariate distribution (median [IQR]):",
        *[f"  {c:20s} {df[c].median():.3f} [{df[c].quantile(.25):.3f}-{df[c].quantile(.75):.3f}]" for c in cov],
        "",
        "Spearman correlation with the v1 SMR among municipalities with >= 50 deaths:",
        *[f"  {k:20s} {v:+.3f}" for k, v in corr.items()],
    ]
    (V2 / "data" / "dataset_report.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
