#!/usr/bin/env python3
"""
vitals.py -- Python port of 05_vitals.sql.

Output microbiology_Vitals.csv: one row per
(anon_id, pat_enc_csn_id_coded, order_proc_id_coded), from the 48h before each
culture order. Per vital, 5 columns: Q25/Q75/median/first/last. Vitals:
heartrate, resprate, temp, sysbp, diasbp. (weight is filtered but never
aggregated in the SQL, so it is omitted here too.)

Run: python vitals.py  -- edit the CONFIGURATION block first.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# ============================================================================
# CONFIGURATION
# ============================================================================

RAW_DIR    = "/PATH/TO/raw/edw/2024-12-31"
COHORT_CSV = "/PATH/TO/microbiology_cultures_cohort.csv"
OUTPUT_CSV = "microbiology_Vitals.csv"

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

# ADAPT: one CSV per vital. hr/resp/temp share a single-value schema; bp holds
# bpSystolic + bpDiastolic columns. Fix any filename that differs.
HRTEMP_FILES: Dict[str, str] = {
    "heartrate": os.path.join(RAW_DIR, "heartrate.csv"),
    "resprate":  os.path.join(RAW_DIR, "resprate.csv"),
    "temp":      os.path.join(RAW_DIR, "tempt.csv"),
}
BP_FILE = os.path.join(RAW_DIR, "bp.csv")

# Lookback window in hours (matches SQL).
LOOKBACK_HOURS = 48


# ============================================================================
# SQL -> MGB SCHEMA MAPPING
#
# Stanford BigQuery                   MGB EDW CSV
# ---------------------------------   ----------------------------------------
# c.anon_id                           patientID (raw) / anon_id (deid)
# c.pat_enc_csn_id_coded              patientEncounterID (raw) / ..._coded (deid)
# c.order_proc_id_coded               orderProcedureID
# c.order_time_jittered_utc           cultureCollDTS
# v.recorded_time_jittered_utc        recordedDTS
# v.row_disp_name                     implicit from FILE NAME (one vital / file)
# v.meas_value (numeric string)       resultN (hr/resp/temp)
# v.meas_value ("120/80" for BP)      bpSystolic + bpDiastolic (TWO columns)
#
# JOIN: the SQL joins vitals on v.anon_id = c.anon_id ONLY (not encounter), so
# a vital within 48h attaches to a culture regardless of encounter. We keep
# this; the output row's encounter is the cohort's, not the vital's.
#
# ADAPT checklist -- verify against your MGB export:
# (A) One file per vital, so no row_disp_name filter is needed. If a file mixes
#     types, add a FLOWSHEET_ID_FILTERS entry.
# (B) BP: SQL stores "SYS/DIA"; MGB stores numeric bpSystolic + bpDiastolic.
#     read_bp() handles it -- only touch it if the column names differ.
# (C) Temp unit: SQL's `>= 90` implies Fahrenheit. Confirm MGB temp is in F;
#     if Celsius, convert upstream or change TEMP_MIN (~32).
#
# Deliberate matches to the SQL (do NOT "fix"):
# * hr/resp/temp/bp pulled patient-only (no encounter constraint).
# * weight filtered but not aggregated => no weight columns emitted.
# ============================================================================

# (A) ADAPT: allowed flowsheet IDs if a file mixes measurement types.
# e.g. FLOWSHEET_ID_FILTERS["heartrate"] = {"301070"}
FLOWSHEET_ID_FILTERS: Dict[str, set] = {
}

# (C) Value-range filters (verbatim from the SQL WHERE clause).
HEARTRATE_MIN = 30
RESPRATE_MIN  = 4
RESPRATE_MAX  = 60
TEMP_MIN      = 90     # Fahrenheit; see checklist item (C)
SYSBP_MIN     = 40
DIASBP_MIN    = 30

# ADAPT: column names in hr/resp/temp CSVs.
V_ID_COL       = "patientID"
V_TIME_COL     = "recordedDTS"
V_VALUE_N_COL  = "resultN"
V_VALUE_TXT_COL = "measureTXT"   # optional fallback
V_FMID_COL     = "flowsheetMeasureID"

# ADAPT: column names in bp.csv (two numeric columns, not "SYS/DIA" strings).
BP_ID_COL   = "patientID"
BP_TIME_COL = "recordedDTS"
BP_SYS_COL  = "bpSystolic"
BP_DIA_COL  = "bpDiastolic"


# ============================================================================
# I/O
# ============================================================================

def read_cohort() -> pd.DataFrame:
    df = pd.read_csv(COHORT_CSV, low_memory=False)
    required = [COHORT_PATIENT_COL, COHORT_ENCOUNTER_COL,
                COHORT_ORDER_COL, COHORT_TIME_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Cohort missing columns: {missing}. "
            "Fix COHORT_*_COL constants at the top of the script."
        )
    out = df[required].rename(columns={
        COHORT_PATIENT_COL:   "anon_id",
        COHORT_ENCOUNTER_COL: "pat_enc_csn_id_coded",
        COHORT_ORDER_COL:     "order_proc_id_coded",
        COHORT_TIME_COL:      "order_time_jittered_utc",
    })
    out["order_time_jittered_utc"] = pd.to_datetime(
        out["order_time_jittered_utc"], errors="coerce", utc=False
    )
    out["anon_id"]              = out["anon_id"].astype(str)
    out["pat_enc_csn_id_coded"] = out["pat_enc_csn_id_coded"].astype(str)
    out["order_proc_id_coded"]  = out["order_proc_id_coded"].astype(str)
    return out


def _optional_numeric_fallback(primary: pd.Series, txt: Optional[pd.Series]) -> pd.Series:
    v = pd.to_numeric(primary, errors="coerce")
    if txt is None:
        return v
    fb = pd.to_numeric(txt, errors="coerce")
    return v.where(v.notna(), fb)


def read_hrtemp(label: str, path: str,
                lo: Optional[float] = None,
                hi: Optional[float] = None) -> Optional[pd.DataFrame]:
    if not os.path.exists(path):
        print(f"  [skip] {label}: file not found at {path}", file=sys.stderr)
        return None
    df = pd.read_csv(path, low_memory=False)
    required = [V_ID_COL, V_TIME_COL, V_VALUE_N_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")

    ids = FLOWSHEET_ID_FILTERS.get(label)
    if ids and V_FMID_COL in df.columns:
        df = df[df[V_FMID_COL].astype(str).isin(ids)]

    txt = df[V_VALUE_TXT_COL] if V_VALUE_TXT_COL in df.columns else None
    # SQL filters on the RAW cast, then rounds only in projection. We match:
    # filter on raw, round only what we emit.
    val_raw = _optional_numeric_fallback(df[V_VALUE_N_COL], txt)

    keep = val_raw.notna()
    if lo is not None:
        keep &= val_raw >= lo
    if hi is not None:
        keep &= val_raw <= hi
    df = df[keep].copy()
    if df.empty:
        return None

    return pd.DataFrame({
        "anon_id":   df[V_ID_COL].astype(str),
        "vitaltime": pd.to_datetime(df[V_TIME_COL], errors="coerce"),
        "value":     val_raw[keep].round(2).values,
    }).dropna(subset=["vitaltime"])


def read_bp(path: str) -> Optional[Dict[str, pd.DataFrame]]:
    if not os.path.exists(path):
        print(f"  [skip] bp: file not found at {path}", file=sys.stderr)
        return None
    df = pd.read_csv(path, low_memory=False)
    required = [BP_ID_COL, BP_TIME_COL, BP_SYS_COL, BP_DIA_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")

    ids = FLOWSHEET_ID_FILTERS.get("bp")
    if ids and V_FMID_COL in df.columns:
        df = df[df[V_FMID_COL].astype(str).isin(ids)]

    # Filter on RAW cast (SQL semantics), round on emit.
    sys_raw = pd.to_numeric(df[BP_SYS_COL], errors="coerce")
    dia_raw = pd.to_numeric(df[BP_DIA_COL], errors="coerce")
    t       = pd.to_datetime(df[BP_TIME_COL], errors="coerce")
    pid     = df[BP_ID_COL].astype(str)

    m_sys = sys_raw.notna() & (sys_raw >= SYSBP_MIN) & t.notna()
    m_dia = dia_raw.notna() & (dia_raw >= DIASBP_MIN) & t.notna()

    return {
        "sysbp":  pd.DataFrame({"anon_id": pid[m_sys], "vitaltime": t[m_sys],
                                "value": sys_raw[m_sys].round(2)}),
        "diasbp": pd.DataFrame({"anon_id": pid[m_dia], "vitaltime": t[m_dia],
                                "value": dia_raw[m_dia].round(2)}),
    }


# ============================================================================
# 48h WINDOWED AGGREGATION -- per vital stream, join to cohort on anon_id
# alone (patient-only, per SQL), keep TIMESTAMP_DIFF(vital, cohort, HOUR) in
# [-48, 0], then aggregate per (anon_id, pat_enc_csn_id_coded, order_proc_id_coded).
# ============================================================================

def _q(x: pd.Series, p: float) -> float:
    v = x.dropna()
    if v.empty:
        return np.nan
    return float(np.quantile(v, p))


def _agg_one_vital(events: pd.DataFrame, cohort: pd.DataFrame, label: str) -> pd.DataFrame:
    if events is None or events.empty:
        return pd.DataFrame()

    j = cohort.merge(events, on="anon_id", how="inner")
    if j.empty:
        return pd.DataFrame()

    # TIMESTAMP_DIFF(..., HOUR) truncates toward zero, so a vital 48.5h before
    # the culture yields -48 and is kept by BETWEEN -48 AND 0. np.trunc matches.
    dt_sec = (j["vitaltime"] - j["order_time_jittered_utc"]).dt.total_seconds()
    j["dt_hr"] = np.trunc(dt_sec / 3600.0).astype("Int64")
    j = j[(j["dt_hr"] >= -LOOKBACK_HOURS) & (j["dt_hr"] <= 0)]
    if j.empty:
        return pd.DataFrame()

    j = j.sort_values(
        ["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded", "vitaltime"],
        kind="mergesort",   # stable: matches R data.table setorder tie-breaking
    )

    def _summarize(g: pd.DataFrame) -> pd.Series:
        v = g["value"]
        return pd.Series({
            f"Q25_{label}":    round(_q(v, 0.25), 2),
            f"Q75_{label}":    round(_q(v, 0.75), 2),
            f"median_{label}": round(_q(v, 0.50), 2),
            f"first_{label}":  round(float(v.iloc[0]),  2),
            f"last_{label}":   round(float(v.iloc[-1]), 2),
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


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print(f"Reading cohort: {COHORT_CSV}", file=sys.stderr)
    cohort = read_cohort()
    print(f"Cohort rows: {len(cohort):,}", file=sys.stderr)

    streams: Dict[str, pd.DataFrame] = {}

    print("\nReading hr / resp / temp ...", file=sys.stderr)
    streams["heartrate"] = read_hrtemp("heartrate", HRTEMP_FILES["heartrate"],
                                       lo=HEARTRATE_MIN)
    streams["resprate"]  = read_hrtemp("resprate",  HRTEMP_FILES["resprate"],
                                       lo=RESPRATE_MIN, hi=RESPRATE_MAX)
    streams["temp"]      = read_hrtemp("temp",      HRTEMP_FILES["temp"],
                                       lo=TEMP_MIN)

    print("\nReading bp ...", file=sys.stderr)
    bp = read_bp(BP_FILE)
    if bp is not None:
        streams["sysbp"]  = bp["sysbp"]
        streams["diasbp"] = bp["diasbp"]

    print("\nExtracting 48h features ...", file=sys.stderr)
    parts: List[pd.DataFrame] = []
    for label in ("heartrate", "resprate", "temp", "sysbp", "diasbp"):
        s = streams.get(label)
        if s is None or s.empty:
            print(f"  [skip] {label}: empty stream", file=sys.stderr)
            continue
        print(f"  + {label}", file=sys.stderr, flush=True)
        agg = _agg_one_vital(s, cohort, label)
        if not agg.empty:
            parts.append(agg)

    if not parts:
        print("No output rows.", file=sys.stderr); return

    out = parts[0]
    for p in parts[1:]:
        out = out.merge(
            p,
            on=["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded"],
            how="outer",
        )

    # Column order: keys first, then per-vital metrics grouped by vital.
    keys = ["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded"]
    ordered = list(keys)
    for label in ("heartrate", "resprate", "temp", "sysbp", "diasbp"):
        for stat in ("Q25", "Q75", "median", "first", "last"):
            col = f"{stat}_{label}"
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
