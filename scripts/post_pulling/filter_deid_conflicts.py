"""
filter_deid_conflicts.py
------------------------
Filters 'different_patients' conflict rows from an already-deidentified data
pull, where the raw CONFLICT_...csv file is not directly applicable (its
order_proc_id_coded column contains raw IDs like 'BWFH_F...' but the pull uses
deidentified integer IDs).

Detection: any order_proc_id_coded that appears with more than one distinct
anon_id in the cohort file is a 'different_patients' conflict. Those orders
are removed from every CSV in the input folder that has an
order_proc_id_coded column.

Usage:
    python filter_deid_conflicts.py --input <folder> --output <folder>
        [--cohort microbiology_cultures_cohort.csv]
"""

import argparse
import os
import sys

import pandas as pd

DEFAULT_COHORT_NAME = "microbiology_cultures_cohort.csv"


def find_conflict_orders(cohort_path: str) -> set:
    """Return order_proc_id_coded values that map to >1 anon_id (different_patients)."""
    df = pd.read_csv(cohort_path, usecols=["order_proc_id_coded", "anon_id"])
    counts = df.drop_duplicates().groupby("order_proc_id_coded")["anon_id"].nunique()
    bad = set(counts[counts > 1].index)
    print(f"Detected {len(bad):,} 'different_patients' order_proc_id_coded values.")
    return bad


def filter_file(path: str, bad_ids: set, out_path: str) -> dict:
    df = pd.read_csv(path, low_memory=False)
    n_before = len(df)

    if "order_proc_id_coded" not in df.columns:
        df.to_csv(out_path, index=False)
        return {"rows_before": n_before, "rows_removed": 0, "rows_after": n_before}

    mask = df["order_proc_id_coded"].isin(bad_ids)
    df_filtered = df[~mask]
    df_filtered.to_csv(out_path, index=False)
    return {
        "rows_before": n_before,
        "rows_removed": int(mask.sum()),
        "rows_after": len(df_filtered),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input folder with deidentified CSVs.")
    parser.add_argument("--output", required=True, help="Output folder for filtered CSVs.")
    parser.add_argument("--cohort", default=DEFAULT_COHORT_NAME,
                        help="Cohort file name (used to detect conflict orders).")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    cohort_path = os.path.join(input_dir, args.cohort)

    if not os.path.isdir(input_dir):
        print(f"ERROR: Input folder not found: {input_dir}")
        sys.exit(1)
    if not os.path.isfile(cohort_path):
        print(f"ERROR: Cohort file not found: {cohort_path}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    bad_ids = find_conflict_orders(cohort_path)

    csv_files = sorted(
        f for f in os.listdir(input_dir)
        if f.lower().endswith(".csv") and not f.startswith("CONFLICT_")
    )

    if not csv_files:
        print(f"No CSV files found in: {input_dir}")
        sys.exit(1)

    print(f"\nProcessing {len(csv_files)} file(s) from: {input_dir}")
    print(f"Writing filtered files to:              {output_dir}\n")

    total_removed = 0
    for fname in csv_files:
        in_path = os.path.join(input_dir, fname)
        out_path = os.path.join(output_dir, fname)
        summary = filter_file(in_path, bad_ids, out_path)
        total_removed += summary["rows_removed"]
        tag = f"removed {summary['rows_removed']:>6,} row(s)" if summary["rows_removed"] else "no rows removed"
        print(f"  {fname:<65s}  {summary['rows_before']:>10,} → {summary['rows_after']:>10,}  |  {tag}")

    print(f"\nDone. Total rows removed across all files: {total_removed:,}")


if __name__ == "__main__":
    main()
