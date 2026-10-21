"""Maternal composition and place of birth, per municipality of residence, 2014-2024.

Reads the eleven yearly SINASC files (one row per live birth) and aggregates, per six-digit
municipality of residence, the counts the model needs:

    composition of the mothers   age < 20, age >= 35, schooling <= 7 years, non-white
    place of birth               hospital birth, birth outside the municipality of
                                 residence, birth in a municipality that held NICU beds
                                 in that year (CNES)
    intermediates (descriptive)  fewer than 7 antenatal visits, preterm (< 37 weeks),
                                 birth weight < 2,500 g

Every share is reported with its own known-denominator so that missing codes (9, blank)
do not count as absence; the known-denominators also give the share of births with missing
schooling or race/colour, used as registration indicators.

Input (individual level, not distributed): SINASC_<year>.csv, semicolon-separated, from DATASUS.
Set NNM1_SINASC to the folder that holds them.

Outputs (data/):
    maternal_composition_muni.csv        pooled 2014-2024, one row per municipality
    maternal_composition_muni_year.csv   the same by year
    maternal_composition_report.txt      totals and reconciliation against the births table

Run (global Python, needs polars):
    python src/01_maternal_composition.py
"""

from pathlib import Path
import time

import polars as pl

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
A1 = V2.parent
SINASC = Path(__import__("os").environ.get("NNM1_SINASC", A1 / "data" / "raw" / "SINASC"))
YEARS = range(2014, 2025)

COLS = ["CODMUNRES", "CODMUNNASC", "LOCNASC", "IDADEMAE", "ESCMAE", "CONSULTAS",
        "GESTACAO", "PESO", "RACACORMAE"]


def nicu_municipalities_by_year():
    """Six-digit codes of municipalities with at least one neonatal ICU bed, per year."""
    cn = pl.read_csv(A1 / "data" / "processed" / "cnes_muni_ano.csv",
                     schema_overrides={"CODMUNRES": pl.Utf8})
    cn = cn.with_columns(pl.col("CODMUNRES").str.slice(0, 6).alias("cod6"))
    out = {}
    for y in YEARS:
        s = cn.filter((pl.col("ANO") == y) & (pl.col("uti_neo_leitos") > 0))
        out[y] = set(s["cod6"].to_list())
    return out


def aggregate_year(year, nicu_set):
    path = SINASC / f"SINASC_{year}.csv"
    lf = pl.scan_csv(path, separator=";", infer_schema_length=0, ignore_errors=True)
    lf = lf.rename({c: c.upper() for c in lf.collect_schema().names()})  # some years use lower-case headers
    lf = lf.select(COLS).with_columns(pl.col(c).cast(pl.Utf8) for c in COLS)
    age = pl.col("IDADEMAE").cast(pl.Int32, strict=False)
    peso = pl.col("PESO").cast(pl.Int32, strict=False)
    esc = pl.col("ESCMAE")
    cons = pl.col("CONSULTAS")
    gest = pl.col("GESTACAO")
    race = pl.col("RACACORMAE")
    loc = pl.col("LOCNASC")
    res = pl.col("CODMUNRES").str.slice(0, 6)
    nasc = pl.col("CODMUNNASC").str.slice(0, 6)

    lf = lf.with_columns([
        res.alias("cod6"),
        pl.lit(1).alias("births"),
        age.is_between(10, 60).cast(pl.Int32).alias("age_known"),
        (age < 20).fill_null(False).cast(pl.Int32).alias("age_lt20"),
        (age >= 35).fill_null(False).cast(pl.Int32).alias("age_ge35"),
        esc.is_in(["1", "2", "3", "4", "5"]).cast(pl.Int32).alias("esc_known"),
        esc.is_in(["1", "2", "3"]).cast(pl.Int32).alias("esc_le7"),
        esc.is_in(["1", "2"]).cast(pl.Int32).alias("esc_le3"),
        cons.is_in(["1", "2", "3", "4"]).cast(pl.Int32).alias("cons_known"),
        cons.is_in(["1", "2", "3"]).cast(pl.Int32).alias("cons_lt7"),
        gest.is_in(["1", "2", "3", "4", "5", "6"]).cast(pl.Int32).alias("gest_known"),
        gest.is_in(["1", "2", "3", "4"]).cast(pl.Int32).alias("preterm"),
        peso.is_between(100, 8000).cast(pl.Int32).alias("peso_known"),
        ((peso < 2500) & peso.is_between(100, 8000)).fill_null(False).cast(pl.Int32).alias("lbw"),
        race.is_in(["1", "2", "3", "4", "5"]).cast(pl.Int32).alias("race_known"),
        race.is_in(["2", "4", "5"]).cast(pl.Int32).alias("race_nonwhite"),
        loc.is_in(["1", "2", "3", "4", "5"]).cast(pl.Int32).alias("loc_known"),
        (loc == "1").cast(pl.Int32).alias("loc_hospital"),
        (nasc.is_not_null() & (nasc != res)).cast(pl.Int32).alias("born_outside"),
        nasc.is_not_null().cast(pl.Int32).alias("nasc_known"),
        nasc.is_in(sorted(nicu_set)).cast(pl.Int32).alias("born_in_nicu_muni"),
    ])
    sums = [c for c in ["births", "age_known", "age_lt20", "age_ge35", "esc_known", "esc_le7", "esc_le3",
                        "cons_known", "cons_lt7", "gest_known", "preterm", "peso_known", "lbw",
                        "race_known", "race_nonwhite", "loc_known", "loc_hospital", "born_outside",
                        "nasc_known", "born_in_nicu_muni"]]
    agg = lf.group_by("cod6").agg([pl.col(c).sum() for c in sums]).with_columns(pl.lit(year).alias("year"))
    return agg.collect()


def main():
    t0 = time.time()
    nicu = nicu_municipalities_by_year()
    frames = []
    for y in YEARS:
        t = time.time()
        f = aggregate_year(y, nicu[y])
        frames.append(f)
        print(f"{y}: {f['births'].sum():,} births, {f.height:,} municipalities, {time.time() - t:.0f}s", flush=True)
    by_year = pl.concat(frames).sort(["cod6", "year"])
    by_year.write_csv(V2 / "data" / "maternal_composition_muni_year.csv")
    pooled = by_year.drop("year").group_by("cod6").agg(pl.all().sum()).sort("cod6")
    pooled = pooled.rename({"cod6": "CODMUNRES"})
    pooled.write_csv(V2 / "data" / "maternal_composition_muni.csv")

    # Reconciliation against the v1 denominators.
    v1_births = pl.read_csv(A1 / "data" / "processed" / "nascidos_muni_ano.csv",
                            schema_overrides={"CODMUNRES": pl.Utf8})
    v1_births = v1_births.with_columns(pl.col("CODMUNRES").str.slice(0, 6).alias("cod6"))
    v1_tot = v1_births["nascimentos"].sum()
    v1_m = v1_births.group_by("cod6").agg(pl.col("nascimentos").sum())
    cmp = pooled.join(v1_m, left_on="CODMUNRES", right_on="cod6", how="full")
    diff = (cmp["births"].fill_null(0) - cmp["nascimentos"].fill_null(0))
    tot = pooled.sum()
    lines = [
        "Maternal composition, pooled 2014-2024",
        "=" * 60,
        f"births in SINASC files             : {tot['births'][0]:,}",
        f"births in v1 nascidos_muni_ano.csv : {v1_tot:,}",
        f"municipality-level |difference|    : max {diff.abs().max():,}  sum {diff.abs().sum():,}",
        f"municipalities (pooled)            : {pooled.height:,}",
        "",
        "national shares (known denominators):",
        f"  mother < 20 y                    : {tot['age_lt20'][0] / tot['age_known'][0]:.3f}",
        f"  mother >= 35 y                   : {tot['age_ge35'][0] / tot['age_known'][0]:.3f}",
        f"  schooling <= 7 y                 : {tot['esc_le7'][0] / tot['esc_known'][0]:.3f}",
        f"  schooling <= 3 y                 : {tot['esc_le3'][0] / tot['esc_known'][0]:.3f}",
        f"  non-white mother                 : {tot['race_nonwhite'][0] / tot['race_known'][0]:.3f}",
        f"  < 7 antenatal visits             : {tot['cons_lt7'][0] / tot['cons_known'][0]:.3f}",
        f"  preterm (< 37 w)                 : {tot['preterm'][0] / tot['gest_known'][0]:.3f}",
        f"  birth weight < 2,500 g           : {tot['lbw'][0] / tot['peso_known'][0]:.3f}",
        f"  hospital birth                   : {tot['loc_hospital'][0] / tot['loc_known'][0]:.3f}",
        f"  born outside municipality        : {tot['born_outside'][0] / tot['nasc_known'][0]:.3f}",
        f"  born in a NICU municipality      : {tot['born_in_nicu_muni'][0] / tot['nasc_known'][0]:.3f}",
        "",
        f"missingness: age {1 - tot['age_known'][0] / tot['births'][0]:.4f}, schooling {1 - tot['esc_known'][0] / tot['births'][0]:.4f}, "
        f"race {1 - tot['race_known'][0] / tot['births'][0]:.4f}, visits {1 - tot['cons_known'][0] / tot['births'][0]:.4f}, "
        f"gestation {1 - tot['gest_known'][0] / tot['births'][0]:.4f}, weight {1 - tot['peso_known'][0] / tot['births'][0]:.4f}",
        f"elapsed {time.time() - t0:.0f}s",
    ]
    (V2 / "data" / "maternal_composition_report.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
