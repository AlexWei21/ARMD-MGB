#!/usr/bin/env python3
"""
implied_susceptibility.py
--------------------------
Translates the implied susceptibility SQL logic from Implied-susceptibility.sql
into a pandas-based processing script.

For each organism-antibiotic pair where susceptibility was not directly measured,
an implied susceptibility is inferred from related antibiotics tested on the same
culture order (same order_proc_id_coded + organism).

Organisms covered:
    Acinetobacter, Citrobacter, Enterobacter, Escherichia coli,
    Klebsiella, Proteus, Pseudomonas aeruginosa, Serratia

Input CSV must contain columns:
    anon_id, pat_enc_csn_id_coded, order_proc_id_coded, organism, antibiotic, susceptibility

Usage:
    python implied_susceptibility.py \\
        --input  data/repulled_data/microbiology_cohort.csv \\
        --output data/repulled_data/microbiology_implied_susceptibility.csv
"""

import argparse
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Rule definitions
# Each entry describes one block of implied-susceptibility logic for a
# given organism pattern (SQL LIKE expression translated to a substring).
#
# Rule types
# ----------
# "unconditional_resistant"
#     antibiotics listed are always set to Resistant when susceptibility is NULL.
#
# "conditional_equiv"
#     target_abx gets the value (Susceptible or Resistant) of source_abx when
#     target_abx susceptibility is NULL.
#     When multiple source_abx are listed, each is evaluated in order; the first
#     match (non-null result) wins.
#
# "conditional_susceptible_if_any"
#     target_abx is set to Susceptible if ANY of source_abx is Susceptible for
#     the same order; otherwise no inference.
#
# "conditional_susceptible_if_any_else_resistant"
#     Like conditional_susceptible_if_any but falls back to Resistant if no
#     source is Susceptible.
# ---------------------------------------------------------------------------

RULES = [
    # =====================================================================
    # ACINETOBACTER
    # =====================================================================
    {
        "organism_pattern": "ACINETOBACTER",
        "type": "unconditional_resistant",
        "antibiotics": ["aztreonam", "cefazolin", "minocycline", "tetracycline"],
    },
    {
        "organism_pattern": "ACINETOBACTER",
        "type": "conditional_equiv",
        "target_abx": "doripenem",
        "source_abx": ["meropenem"],
    },
    {
        "organism_pattern": "ACINETOBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "meropenem",
        "source_abx": ["imipenem"],
    },

    # =====================================================================
    # CITROBACTER
    # =====================================================================
    {
        "organism_pattern": "CITROBACTER",
        "type": "unconditional_resistant",
        "antibiotics": ["ampicillin", "cefazolin", "cefotetan", "cefoxitin"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "cefepime",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "ceftazidime",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_equiv",
        "target_abx": "ceftriaxone",
        "source_abx": ["cefotaxime"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_equiv",
        "target_abx": "cefotaxime",
        "source_abx": ["ceftriaxone"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_equiv",
        "target_abx": "doripenem",
        "source_abx": ["meropenem"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "ertapenem",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "imipenem",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "CITROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "meropenem",
        "source_abx": ["imipenem", "ceftriaxone", "cefotaxime"],
    },

    # =====================================================================
    # ENTEROBACTER
    # =====================================================================
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "unconditional_resistant",
        "antibiotics": ["ampicillin", "cefazolin"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "cefepime",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "ceftazidime",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_equiv",
        "target_abx": "ceftriaxone",
        "source_abx": ["cefotaxime"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_equiv",
        "target_abx": "cefotaxime",
        "source_abx": ["ceftriaxone"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_susceptible_if_any_else_resistant",
        "target_abx": "doxycycline",
        "source_abx": ["tetracycline"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_equiv",
        "target_abx": "doripenem",
        "source_abx": ["meropenem"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "ertapenem",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "imipenem",
        "source_abx": ["ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "ENTEROBACTER",
        "type": "conditional_susceptible_if_any",
        "target_abx": "meropenem",
        "source_abx": ["imipenem", "ceftriaxone", "cefotaxime"],
    },

    # =====================================================================
    # ESCHERICHIA COLI  (pattern matches ESCHERICHIA.*COLI)
    # =====================================================================
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_susceptible_if_any",
        "target_abx": "cefepime",
        "source_abx": ["cefazolin", "ceftriaxone", "cefotaxime"],
    },
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_susceptible_if_any",
        "target_abx": "ceftazidime",
        "source_abx": ["cefazolin", "ceftriaxone", "cefotaxime"],
    },
    # ceftriaxone: first equiv cefotaxime, then fallback to cefazolin-S
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_equiv_then_susceptible_if_any",
        "target_abx": "ceftriaxone",
        "equiv_source": "cefotaxime",
        "fallback_source_abx": ["cefazolin"],
    },
    # cefotaxime: first equiv ceftriaxone, then fallback to cefazolin-S
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_equiv_then_susceptible_if_any",
        "target_abx": "cefotaxime",
        "equiv_source": "ceftriaxone",
        "fallback_source_abx": ["cefazolin"],
    },
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_susceptible_if_any",
        "target_abx": "cefotetan",
        "source_abx": ["cefazolin"],
    },
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_susceptible_if_any",
        "target_abx": "cefoxitin",
        "source_abx": ["cefazolin"],
    },
    # cefpodoxime: equiv ceftriaxone, then equiv cefotaxime, then cefazolin-S
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_equiv_chain_then_susceptible_if_any",
        "target_abx": "cefpodoxime",
        "equiv_chain": ["ceftriaxone", "cefotaxime"],
        "fallback_source_abx": ["cefazolin"],
    },
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_susceptible_if_any",
        "target_abx": "cefuroxime",
        "source_abx": ["cefazolin"],
    },
    {
        "organism_pattern": "ESCHERICHIA COLI",
        "type": "conditional_susceptible_if_any_else_resistant",
        "target_abx": "doxycycline",
        "source_abx": ["tetracycline"],
    },

    # =====================================================================
    # KLEBSIELLA
    # =====================================================================
    {
        "organism_pattern": "KLEBSIELLA",
        "type": "unconditional_resistant",
        "antibiotics": ["ampicillin"],
    },

    # =====================================================================
    # PROTEUS
    # =====================================================================
    {
        "organism_pattern": "PROTEUS",
        "type": "unconditional_resistant",
        "antibiotics": ["tetracycline", "tigecycline", "colistin"],
    },

    # =====================================================================
    # PSEUDOMONAS
    # =====================================================================
    {
        "organism_pattern": "PSEUDOMONAS",
        "type": "unconditional_resistant",
        "antibiotics": [
            "ampicillin", "ceftriaxone", "cefazolin", "ertapenem",
            "tetracycline", "tigecycline",
        ],
    },

    # =====================================================================
    # SERRATIA
    # =====================================================================
    {
        "organism_pattern": "SERRATIA",
        "type": "unconditional_resistant",
        "antibiotics": ["ampicillin", "cefazolin", "tetracycline"],
    },
]


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

ID_COLS = ["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded", "organism"]


def organism_mask(df: pd.DataFrame, pattern: str) -> pd.Series:
    """Return a boolean mask for rows where organism matches the pattern."""
    pat = pattern.upper()
    org = df["organism"].str.upper()
    if " " in pat:
        # Multi-word: build a single regex (e.g. ESCHERICHIA.*COLI)
        parts = pat.split()
        regex = ".*".join(parts)
        return org.str.contains(regex, na=False, regex=True)
    return org.str.contains(pat, na=False, regex=False)


def ensure_rows(df: pd.DataFrame, mask: pd.Series, antibiotics: list[str]) -> pd.DataFrame:
    """
    For every (anon_id, pat_enc_csn_id_coded, order_proc_id_coded, organism)
    group matched by mask, ensure a row exists for each antibiotic in the list.
    Missing rows are appended with NULL susceptibility / implied_susceptibility.
    """
    orders = df.loc[mask, ID_COLS].drop_duplicates()
    new_rows = []
    for abx in antibiotics:
        existing = df.loc[
            mask & (df["antibiotic"].str.lower() == abx),
            ID_COLS,
        ].drop_duplicates()
        # rows in orders but not in existing
        missing = orders.merge(existing, on=ID_COLS, how="left", indicator=True)
        missing = missing.loc[missing["_merge"] == "left_only", ID_COLS]
        if not missing.empty:
            missing = missing.copy()
            missing["antibiotic"] = abx
            missing["susceptibility"] = None
            missing["implied_susceptibility"] = None
            new_rows.append(missing)
    if new_rows:
        df = pd.concat([df] + new_rows, ignore_index=True)
    return df


def get_order_set(df: pd.DataFrame, mask: pd.Series, abx_list: list[str],
                  result: str) -> set:
    """
    Return the set of order_proc_id_coded values where organism matches mask
    AND antibiotic (lower) is in abx_list AND (susceptibility OR
    implied_susceptibility) == result.
    """
    abx_mask = df["antibiotic"].str.lower().isin(abx_list)
    result_mask = (
        (df["susceptibility"] == result) |
        (df["implied_susceptibility"] == result)
    )
    return set(df.loc[mask & abx_mask & result_mask, "order_proc_id_coded"])


def get_order_set_suscept_only(df: pd.DataFrame, mask: pd.Series,
                               abx_list: list[str], result: str) -> set:
    """
    Same as get_order_set but checks only susceptibility (not implied)
    – used by rules that mirror the SQL WHERE susceptibility = 'Susceptible'.
    """
    abx_mask = df["antibiotic"].str.lower().isin(abx_list)
    result_mask = df["susceptibility"] == result
    return set(df.loc[mask & abx_mask & result_mask, "order_proc_id_coded"])


def apply_update(df: pd.DataFrame, row_mask: pd.Series,
                 order_ids: set, value: str) -> pd.DataFrame:
    """Set implied_susceptibility = value where row_mask AND order is in order_ids."""
    target = (
        row_mask
        & df["order_proc_id_coded"].isin(order_ids)
        & df["susceptibility"].isna()
        & df["implied_susceptibility"].isna()
    )
    df.loc[target, "implied_susceptibility"] = value
    return df


# ---------------------------------------------------------------------------
# Rule applicators
# Note: all apply_* functions accept the organism_pattern string and recompute
# the mask after ensure_rows, because pd.concat resets the DataFrame index and
# the old boolean mask would be misaligned.
# ---------------------------------------------------------------------------

def apply_unconditional_resistant(df, pattern, antibiotics):
    for abx in antibiotics:
        omask = organism_mask(df, pattern)
        df = ensure_rows(df, omask, [abx])
        omask = organism_mask(df, pattern)  # recompute after concat
        abx_mask = df["antibiotic"].str.lower() == abx
        target = omask & abx_mask & df["susceptibility"].isna()
        df.loc[target, "implied_susceptibility"] = df.loc[target, "implied_susceptibility"].fillna("Resistant")
    return df


def apply_conditional_equiv(df, pattern, target_abx, source_abx_list):
    """
    target_abx gets the value of the first source_abx (S or R) found for the
    same order_proc_id_coded, when target_abx susceptibility is NULL.
    """
    omask = organism_mask(df, pattern)
    df = ensure_rows(df, omask, [target_abx])
    omask = organism_mask(df, pattern)
    target_row_mask = omask & (df["antibiotic"].str.lower() == target_abx) & df["susceptibility"].isna()

    for result in ["Susceptible", "Resistant"]:
        orders = get_order_set(df, omask, source_abx_list, result)
        df = apply_update(df, target_row_mask, orders, result)
    return df


def apply_conditional_susceptible_if_any(df, pattern, target_abx, source_abx_list,
                                          use_implied_in_source=True):
    """
    target_abx → Susceptible if any source is Susceptible; else no inference.
    """
    omask = organism_mask(df, pattern)
    df = ensure_rows(df, omask, [target_abx])
    omask = organism_mask(df, pattern)
    target_row_mask = omask & (df["antibiotic"].str.lower() == target_abx) & df["susceptibility"].isna()

    if use_implied_in_source:
        sus_orders = get_order_set(df, omask, source_abx_list, "Susceptible")
    else:
        sus_orders = get_order_set_suscept_only(df, omask, source_abx_list, "Susceptible")
    df = apply_update(df, target_row_mask, sus_orders, "Susceptible")
    return df


def apply_conditional_susceptible_if_any_else_resistant(df, pattern, target_abx,
                                                          source_abx_list):
    """
    target_abx → Susceptible if source is Susceptible; else Resistant.
    """
    omask = organism_mask(df, pattern)
    df = ensure_rows(df, omask, [target_abx])
    omask = organism_mask(df, pattern)
    target_row_mask = omask & (df["antibiotic"].str.lower() == target_abx) & df["susceptibility"].isna()

    sus_orders = get_order_set(df, omask, source_abx_list, "Susceptible")
    df = apply_update(df, target_row_mask, sus_orders, "Susceptible")

    # Set remaining NULL to Resistant
    still_null = target_row_mask & df["implied_susceptibility"].isna()
    df.loc[still_null, "implied_susceptibility"] = "Resistant"
    return df


def apply_conditional_equiv_then_susceptible_if_any(df, pattern, target_abx,
                                                     equiv_source, fallback_source_abx):
    """
    target_abx:
      1. Equiv equiv_source (S or R) if available.
      2. Susceptible if any fallback_source_abx is Susceptible.
    """
    omask = organism_mask(df, pattern)
    df = ensure_rows(df, omask, [target_abx])
    omask = organism_mask(df, pattern)
    target_row_mask = omask & (df["antibiotic"].str.lower() == target_abx) & df["susceptibility"].isna()

    for result in ["Susceptible", "Resistant"]:
        orders = get_order_set(df, omask, [equiv_source], result)
        df = apply_update(df, target_row_mask, orders, result)

    # fallback: susceptible if fallback source susceptible
    sus_orders = get_order_set_suscept_only(df, omask, fallback_source_abx, "Susceptible")
    df = apply_update(df, target_row_mask, sus_orders, "Susceptible")
    return df


def apply_conditional_equiv_chain_then_susceptible_if_any(df, pattern, target_abx,
                                                           equiv_chain, fallback_source_abx):
    """
    target_abx:
      1. Equiv equiv_chain[0] if available.
      2. Equiv equiv_chain[1] if available.
      3. Susceptible if any fallback_source_abx is Susceptible.
    """
    omask = organism_mask(df, pattern)
    df = ensure_rows(df, omask, [target_abx])
    omask = organism_mask(df, pattern)
    target_row_mask = omask & (df["antibiotic"].str.lower() == target_abx) & df["susceptibility"].isna()

    for src in equiv_chain:
        for result in ["Susceptible", "Resistant"]:
            orders = get_order_set(df, omask, [src], result)
            df = apply_update(df, target_row_mask, orders, result)

    sus_orders = get_order_set_suscept_only(df, omask, fallback_source_abx, "Susceptible")
    df = apply_update(df, target_row_mask, sus_orders, "Susceptible")
    return df


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def apply_all_rules(df: pd.DataFrame) -> pd.DataFrame:
    """Apply all implied susceptibility rules to a tidy susceptibility DataFrame."""
    if "implied_susceptibility" not in df.columns:
        df["implied_susceptibility"] = None

    for i, rule in enumerate(RULES):
        pattern = rule["organism_pattern"]
        rtype = rule["type"]
        print(f"  Rule {i+1}/{len(RULES)}: {pattern} / {rtype} ...", flush=True)

        if rtype == "unconditional_resistant":
            df = apply_unconditional_resistant(df, pattern, rule["antibiotics"])

        elif rtype == "conditional_equiv":
            df = apply_conditional_equiv(
                df, pattern, rule["target_abx"], rule["source_abx"]
            )

        elif rtype == "conditional_susceptible_if_any":
            df = apply_conditional_susceptible_if_any(
                df, pattern, rule["target_abx"], rule["source_abx"]
            )

        elif rtype == "conditional_susceptible_if_any_else_resistant":
            df = apply_conditional_susceptible_if_any_else_resistant(
                df, pattern, rule["target_abx"], rule["source_abx"]
            )

        elif rtype == "conditional_equiv_then_susceptible_if_any":
            df = apply_conditional_equiv_then_susceptible_if_any(
                df, pattern, rule["target_abx"],
                rule["equiv_source"], rule["fallback_source_abx"]
            )

        elif rtype == "conditional_equiv_chain_then_susceptible_if_any":
            df = apply_conditional_equiv_chain_then_susceptible_if_any(
                df, pattern, rule["target_abx"],
                rule["equiv_chain"], rule["fallback_source_abx"]
            )

        else:
            raise ValueError(f"Unknown rule type: {rtype}")

    return df


def build_implied_susceptibility(df: pd.DataFrame) -> pd.DataFrame:
    """
    Full pipeline:
      1. Normalise antibiotic names to lowercase.
      2. Apply all organism-specific implied susceptibility rules.
      3. Remove rows where both susceptibility and implied_susceptibility are null
         (i.e. retain only informative rows).
      4. Deduplicate.
    Returns the processed DataFrame.
    """
    required = set(ID_COLS + ["antibiotic", "susceptibility"])
    missing_cols = required - set(df.columns)
    if missing_cols:
        raise ValueError(f"Input DataFrame is missing columns: {missing_cols}")

    df = df.copy()
    df["antibiotic"] = df["antibiotic"].str.lower().str.strip()
    df["susceptibility"] = df["susceptibility"].where(df["susceptibility"].notna() & (df["susceptibility"] != ""), None)

    if "implied_susceptibility" not in df.columns:
        df["implied_susceptibility"] = None

    df = apply_all_rules(df)

    # Keep only rows with at least one value
    has_value = df["susceptibility"].notna() | df["implied_susceptibility"].notna()
    df = df[has_value].copy()

    # Deduplicate
    df = df.drop_duplicates(
        subset=["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded",
                "organism", "antibiotic", "susceptibility", "implied_susceptibility"]
    ).sort_values(
        ["anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded", "organism", "antibiotic"]
    ).reset_index(drop=True)

    return df


# ---------------------------------------------------------------------------
# CLI entry-point
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply implied susceptibility rules to a microbiology cohort CSV."
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to the input CSV containing culture susceptibility data. "
             "Required columns: anon_id, pat_enc_csn_id_coded, order_proc_id_coded, "
             "organism, antibiotic, susceptibility."
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Path for the output CSV (will be created or overwritten)."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    print(f"Reading input: {input_path}")
    df = pd.read_csv(input_path, low_memory=False)
    print(f"  Rows loaded: {len(df):,}")

    print("Applying implied susceptibility rules ...")
    result = build_implied_susceptibility(df)
    print(f"  Rows after processing: {len(result):,}")

    print(f"Writing output: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    print("Done.")

    # Summary stats
    total = len(result)
    direct = result["susceptibility"].notna().sum()
    implied_only = (result["susceptibility"].isna() & result["implied_susceptibility"].notna()).sum()
    print(f"\nSummary")
    print(f"  Total rows:         {total:,}")
    print(f"  Direct suscept.:    {direct:,}")
    print(f"  Implied only:       {implied_only:,}")
    print(f"\nImplied susceptibility value counts:")
    print(result["implied_susceptibility"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
