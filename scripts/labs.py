#!/usr/bin/env python3
"""
labs.py -- Python port of 04_labs.sql.

Output microbiology_Labs.csv: one row per
(anon_id, pat_enc_csn_id_coded, order_proc_id_coded, Period_Day),
Period_Day in {14, 30, 180} days. Per lab, 5 columns:
Q25/Q75/median/first/last. Labs: wbc, neutrophils, lymphocytes, hgb, plt,
na, hco3, bun, cr, lactate, procalcitonin. (crp is defined but never
aggregated in the SQL, so it is omitted here too.)

Run: python labs.py  -- edit the CONFIGURATION block first.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


# ============================================================================
# CONFIGURATION
# ============================================================================

RAW_DIR    = "/PATH/TO/raw/edw/2024-12-31"
COHORT_CSV = "/PATH/TO/microbiology_cultures_cohort.csv"
OUTPUT_CSV = "microbiology_Labs.csv"

# ADAPT: deid. Set APPLY_DEID=True to swap raw IDs for the project's deid
# integers before writing. Point MAPS_DIR at the folder holding
# anon_id_map.csv / enc_id_map.csv / order_proc_id_map.csv.
APPLY_DEID = True
MAPS_DIR   = "/PATH/TO/maps"

# ADAPT: cohort column names (your file -> Stanford identifier).
COHORT_PATIENT_COL   = "patientID"           # -> anon_id
COHORT_ENCOUNTER_COL = "patientEncounterID"  # -> pat_enc_csn_id_coded
COHORT_ORDER_COL     = "orderProcedureID"    # -> order_proc_id_coded
COHORT_TIME_COL      = "cultureCollDTS"      # -> order_time_jittered_utc

# ADAPT: one CSV per lab (MGB EDW = one file per variable). Fix any filename
# that differs. crp omitted on purpose (see docstring).
LAB_FILES: Dict[str, str] = {
    "wbc":           os.path.join(RAW_DIR, "wbc.csv"),
    "neutrophils":   os.path.join(RAW_DIR, "neutrophils.csv"),
    "lymphocytes":   os.path.join(RAW_DIR, "lymphocytes.csv"),
    "hgb":           os.path.join(RAW_DIR, "hgb.csv"),
    "plt":           os.path.join(RAW_DIR, "plt.csv"),
    "na":            os.path.join(RAW_DIR, "na.csv"),
    "hco3":          os.path.join(RAW_DIR, "hco3.csv"),
    "bun":           os.path.join(RAW_DIR, "bun.csv"),
    "cr":            os.path.join(RAW_DIR, "cr.csv"),
    "lactate":       os.path.join(RAW_DIR, "lactate.csv"),
    "procalcitonin": os.path.join(RAW_DIR, "procalcitonin.csv"),
}

# Lookback windows in days (matches SQL).
LOOKBACK_DAYS: List[int] = [14, 30, 180]


# ============================================================================
# SQL -> MGB SCHEMA MAPPING
#
# Stanford BigQuery                   MGB EDW CSV
# ---------------------------------   ----------------------------------------
# c.anon_id                           patientID (raw) / anon_id (deid)
# c.pat_enc_csn_id_coded              patientEncounterID (raw) / ..._coded (deid)
# c.order_proc_id_coded               orderProcedureID
# c.order_time_jittered_utc           cultureCollDTS
# lr.order_time_jittered_utc          orderDTS
# lr.base_name / lr.lab_name          implicit from FILE NAME (one lab per file)
# lr.reference_unit                   referenceRangeUnit
# lr.ord_value                        resultN (numeric) / resultTXT (e.g. ">10")
#
# ADAPT checklist -- verify against your MGB export:
# (A) UNIT_MAP is copied verbatim from the SQL (lowercased/trimmed). If your
#     referenceRangeUnit has other spellings/casings (e.g. "10^3/uL"), extend
#     the relevant entry.
# (B) If one CSV mixes lab types, set COMPONENT_NAME_FILTERS[lab] to isolate rows.
# (C) resultN is used first, resultTXT as fallback (USE_RESULTTXT_FALLBACK).
#
# Deliberate matches to the SQL (do NOT "fix"):
# * hgb g/dl -> mg/dl multiplies by 1000 (clinically odd, but matches SQL).
# * hco3 accepts meq/l or mmol/l; lactate accepts base_name 'lac'.
# * neutrophils/lymphocytes require unit '%' (absolute counts dropped).
# ============================================================================

# (A) Unit-string sets; compared lowercased/trimmed.
UNIT_MAP: Dict[str, Dict] = {
    "wbc": {
        "one_to_one":   {"thousand/ul", "k/ul", "10x3/ul", "10*3/ul", "x10e3/ul"},
        "divide_1000":  {"/ul"},   # SQL matches '/uL' case-insensitively
    },
    "neutrophils": {"pct_only": {"%"}},
    "lymphocytes": {"pct_only": {"%"}},
    "hgb": {
        "one_to_one":     {"mg/dl"},
        "multiply_1000":  {"g/dl"},
    },
    "plt": {
        "one_to_one":  {"x10e3/ul", "10x3/ul", "k/ul", "10*3/ul", "thousand/ul"},
        "divide_1000": {"ul"},
    },
    "na":            {"one_to_one": {"mmol/l"}},
    "hco3":          {"one_to_one": {"meq/l", "mmol/l"}},
    "bun":           {"one_to_one": {"mg/dl"}},
    "cr":            {"one_to_one": {"mg/dl"}},
    "lactate":       {"one_to_one": {"mmol/l", "mmole/l"}},
    "procalcitonin": {"one_to_one": {"ng/ml"}},
}

# (B) ADAPT: per-lab componentNM allow-set for files with mixed lab types.
# e.g. COMPONENT_NAME_FILTERS["wbc"] = {"WBC", "WBC AUTO DIFFERENTIAL"}
COMPONENT_NAME_FILTERS: Dict[str, set] = {
}

# (C) True = fall back to to_numeric(resultTXT) when resultN is null.
USE_RESULTTXT_FALLBACK = True

# ADAPT: column names in every MGB lab CSV.
LAB_ID_COL       = "patientID"
LAB_ENC_COL      = "patientEncounterID"
LAB_TIME_COL     = "orderDTS"
LAB_UNIT_COL     = "referenceRangeUnit"
LAB_VALUE_N_COL  = "resultN"
LAB_VALUE_TXT_COL = "resultTXT"
LAB_COMPONENT_COL = "componentNM"


# ============================================================================
# UNIT CONVERTERS -- return the value in the SQL's unified unit, else NaN.
# ============================================================================

def _norm(s: pd.Series) -> pd.Series:
    """Lower + strip (matches SQL LOWER + TRIM)."""
    return s.astype("string").str.strip().str.lower()


def _convert_generic(value: pd.Series, unit: pd.Series, rules: Dict) -> pd.Series:
    """Apply one lab's unit-rule dict (one_to_one/multiply_1000/divide_1000/pct_only)."""
    u = _norm(unit)
    out = pd.Series(np.nan, index=value.index, dtype="float64")
    v = pd.to_numeric(value, errors="coerce")
    if "one_to_one" in rules:
        m = u.isin(rules["one_to_one"])
        out[m] = v[m]
    if "multiply_1000" in rules:
        m = u.isin(rules["multiply_1000"])
        out[m] = v[m] * 1000.0
    if "divide_1000" in rules:
        m = u.isin(rules["divide_1000"])
        out[m] = v[m] / 1000.0
    if "pct_only" in rules:  # same as one_to_one, for '%'
        m = u.isin(rules["pct_only"])
        out[m] = v[m]
    return out


def convert(lab: str, value: pd.Series, unit: pd.Series) -> pd.Series:
    return _convert_generic(value, unit, UNIT_MAP[lab])


# ============================================================================
# I/O
# ============================================================================

def read_cohort() -> pd.DataFrame:
    """Read cohort CSV, rename to Stanford identifiers, coerce time."""
    df = pd.read_csv(COHORT_CSV, low_memory=False)
    required = [COHORT_PATIENT_COL, COHORT_ENCOUNTER_COL,
                COHORT_ORDER_COL, COHORT_TIME_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Cohort missing columns: {missing}. "
            "Fix COHORT_*_COL constants at the top of the script."
        )
    df = df[required].rename(columns={
        COHORT_PATIENT_COL:   "anon_id",
        COHORT_ENCOUNTER_COL: "pat_enc_csn_id_coded",
        COHORT_ORDER_COL:     "order_proc_id_coded",
        COHORT_TIME_COL:      "order_time_jittered_utc",
    })
    df["order_time_jittered_utc"] = pd.to_datetime(
        df["order_time_jittered_utc"], errors="coerce", utc=False
    )
    df["anon_id"]              = df["anon_id"].astype(str)
    df["pat_enc_csn_id_coded"] = df["pat_enc_csn_id_coded"].astype(str)
    df["order_proc_id_coded"]  = df["order_proc_id_coded"].astype(str)
    return df


def read_lab(lab: str, path: str) -> Optional[pd.DataFrame]:
    """Read one lab CSV; return slim frame with unit-converted `value`."""
    if not os.path.exists(path):
        print(f"  [skip] {lab}: file not found at {path}", file=sys.stderr)
        return None
    df = pd.read_csv(path, low_memory=False)
    required = [LAB_ID_COL, LAB_ENC_COL, LAB_TIME_COL, LAB_UNIT_COL, LAB_VALUE_N_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")

    # (B) Optional componentNM filter for files with mixed lab types.
    keep = COMPONENT_NAME_FILTERS.get(lab)
    if keep and LAB_COMPONENT_COL in df.columns:
        df = df[df[LAB_COMPONENT_COL].isin(keep)]

    # (C) resultN first, then optional resultTXT fallback.
    val = pd.to_numeric(df[LAB_VALUE_N_COL], errors="coerce")
    if USE_RESULTTXT_FALLBACK and LAB_VALUE_TXT_COL in df.columns:
        fallback = pd.to_numeric(df[LAB_VALUE_TXT_COL], errors="coerce")
        val = val.where(val.notna(), fallback)

    df = df.assign(value=convert(lab, val, df[LAB_UNIT_COL]))
    df = df[df["value"].notna()].copy()
    if df.empty:
        return None

    return pd.DataFrame({
        "anon_id":              df[LAB_ID_COL].astype(str),
        "pat_enc_csn_id_coded": df[LAB_ENC_COL].astype(str),
        "labtime":              pd.to_datetime(df[LAB_TIME_COL], errors="coerce"),
        "value":                df["value"],
    }).dropna(subset=["labtime"])


# ============================================================================
# WINDOWED AGGREGATION -- per window: join lab stream to cohort on
# (anon_id, pat_enc_csn_id_coded), keep TIMESTAMP_DIFF(lab, cohort, DAY) in
# [-W, 0], then per (patient, enc, order) emit Q25/Q75/median/first/last.
# ============================================================================

def _q(x: pd.Series, p: float) -> float:
    v = x.dropna()
    if v.empty:
        return np.nan
    # np.quantile (linear) approximates SQL APPROX_QUANTILES[OFFSET] to 2dp.
    return float(np.quantile(v, p))


def _agg_one_lab(events: pd.DataFrame, cohort: pd.DataFrame,
                 lookback_days: int, lab: str) -> pd.DataFrame:
    if events is None or events.empty:
        return pd.DataFrame()
    j = cohort.merge(events, on=["anon_id", "pat_enc_csn_id_coded"], how="inner")
    if j.empty:
        return pd.DataFrame()

    # TIMESTAMP_DIFF(..., DAY) truncates toward zero, so a lab 14.5d before the
    # culture yields -14 and is kept by BETWEEN -14 AND 0. np.trunc matches.
    dt_sec = (j["labtime"] - j["order_time_jittered_utc"]).dt.total_seconds()
    j["dt_days"] = np.trunc(dt_sec / 86400.0).astype("Int64")
    j = j[(j["dt_days"] >= -lookback_days) & (j["dt_days"] <= 0)]
    if j.empty:
        return pd.DataFrame()

    j = j.sort_values(
        ["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded", "labtime"],
        kind="mergesort",   # stable: matches R data.table setorder tie-breaking
    )

    def _summarize(g: pd.DataFrame) -> pd.Series:
        v = g["value"]
        return pd.Series({
            f"Q25_{lab}":    round(_q(v, 0.25), 2),
            f"Q75_{lab}":    round(_q(v, 0.75), 2),
            f"median_{lab}": round(_q(v, 0.50), 2),
            f"first_{lab}":  round(float(v.iloc[0]),  2),
            f"last_{lab}":   round(float(v.iloc[-1]), 2),
        })

    agg = (
        j.groupby(
            ["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded"],
            sort=False,
        )
        .apply(_summarize)
        .reset_index()
    )
    return agg


def deidentify_output(df: pd.DataFrame) -> pd.DataFrame:
    """Swap raw IDs for deid surrogates via the three ID maps; warn (don't drop)
    on unmapped rows. No date shift (output has no timestamp columns).
    """
    anon_map = pd.read_csv(os.path.join(MAPS_DIR, "anon_id_map.csv"),
                           dtype={"anon_id": str})
    enc_map  = pd.read_csv(os.path.join(MAPS_DIR, "enc_id_map.csv"),
                           dtype={"pat_enc_csn_id_coded": str})
    ord_map  = pd.read_csv(os.path.join(MAPS_DIR, "order_proc_id_map.csv"),
                           dtype={"order_proc_id_coded": str})

    anon_d = dict(zip(anon_map["anon_id"], anon_map["anon_id_deid"]))
    enc_d  = dict(zip(enc_map["pat_enc_csn_id_coded"],  enc_map["pat_enc_csn_id_deid"]))
    ord_d  = dict(zip(ord_map["order_proc_id_coded"],   ord_map["order_proc_id_deid"]))

    out = df.copy()
    out["anon_id"]              = out["anon_id"].astype(str).map(anon_d)
    out["pat_enc_csn_id_coded"] = out["pat_enc_csn_id_coded"].astype(str).map(enc_d)
    out["order_proc_id_coded"]  = out["order_proc_id_coded"].astype(str).map(ord_d)

    for col in ("anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded"):
        n_missing = int(out[col].isna().sum())
        if n_missing:
            print(f"  WARN: {n_missing:,} rows have no {col} mapping "
                  f"(left as NA)", file=sys.stderr)
    return out


def compute_window(cohort: pd.DataFrame, lookback_days: int,
                   lab_streams: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    print(f"  window: {lookback_days}d", file=sys.stderr)
    parts: List[pd.DataFrame] = []
    for lab, events in lab_streams.items():
        print(f"    + {lab}", file=sys.stderr, flush=True)
        agg = _agg_one_lab(events, cohort, lookback_days, lab)
        if not agg.empty:
            parts.append(agg)
    if not parts:
        return pd.DataFrame()
    out = parts[0]
    for p in parts[1:]:
        out = out.merge(
            p,
            on=["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded"],
            how="outer",
        )
    out["Period_Day"] = lookback_days
    return out


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print(f"Reading cohort: {COHORT_CSV}", file=sys.stderr)
    cohort = read_cohort()
    print(f"Cohort rows: {len(cohort):,}", file=sys.stderr)

    print("\nReading + unit-converting labs...", file=sys.stderr)
    lab_streams: Dict[str, pd.DataFrame] = {}
    for lab, path in LAB_FILES.items():
        s = read_lab(lab, path)
        if s is not None:
            lab_streams[lab] = s
            print(f"  [ok]   {lab:<15s}  {len(s):>10,} rows kept",
                  file=sys.stderr)

    print("\nExtracting windowed features...", file=sys.stderr)
    windows = [compute_window(cohort, w, lab_streams) for w in LOOKBACK_DAYS]
    windows = [w for w in windows if not w.empty]
    if not windows:
        print("No output rows.", file=sys.stderr)
        return

    # UNION ALL + SELECT DISTINCT
    out = pd.concat(windows, ignore_index=True).drop_duplicates()

    # Column order: keys first, then per-lab metrics grouped by lab.
    keys = ["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded", "Period_Day"]
    ordered = list(keys)
    for lab in LAB_FILES.keys():
        for stat in ("Q25", "Q75", "median", "first", "last"):
            col = f"{stat}_{lab}"
            if col in out.columns:
                ordered.append(col)
    out = out[ordered]

    if APPLY_DEID:
        print(f"\nDeidentifying with maps at {MAPS_DIR} ...", file=sys.stderr)
        out = deidentify_output(out)
    else:
        print("\nWARNING: APPLY_DEID = False -- output still contains raw "
              "IDs. Do NOT share externally.", file=sys.stderr)

    out.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out):,} rows -> {OUTPUT_CSV}", file=sys.stderr)


if __name__ == "__main__":
    main()
