"""
filter_different_patients.py
----------------------------
Removes rows from every CSV in a target folder where order_proc_id_coded
appears in the CONFLICT file with mismatch_type == 'different_patients'.

These are genuine source-system data integrity errors (same accession number
issued to two distinct patients at the same hospital) and should be excluded
before any analysis.

Usage:
    python filter_different_patients.py [--conflict PATH] [--input PATH] [--output PATH]

Defaults:
    --conflict  data/new/data/CONFLICT_order_proc_id_coded_to_anon_id_pat_enc_csn_id_coded.csv
    --input     data/new/data/
    --output    data/new/data/filtered/   (created if absent)
"""

import argparse
import os
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Defaults (relative to the repo root)
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

DEFAULT_CONFLICT = os.path.join(
    REPO_ROOT,
    "data", "new", "data",
    "CONFLICT_order_proc_id_coded_to_anon_id_pat_enc_csn_id_coded.csv",
)
DEFAULT_INPUT  = os.path.join(REPO_ROOT, "data", "new", "data")
DEFAULT_OUTPUT = os.path.join(REPO_ROOT, "data", "new", "data", "filtered")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_conflict_ids(conflict_path: str) -> set:
    """Return the set of order_proc_id_coded values flagged as different_patients."""
    df = pd.read_csv(conflict_path)
    bad = df.loc[df["mismatch_type"] == "different_patients", "order_proc_id_coded"]
    ids = set(bad.unique())
    print(f"Loaded {len(ids):,} 'different_patients' order_proc_id_coded values to remove.")
    return ids


def filter_file(path: str, bad_ids: set, out_path: str) -> dict:
    """
    Filter one CSV, return a summary dict.
    Files that don't contain order_proc_id_coded are copied unchanged.
    """
    df = pd.read_csv(path, low_memory=False)
    n_before = len(df)

    if "order_proc_id_coded" not in df.columns:
        df.to_csv(out_path, index=False)
        return {"rows_before": n_before, "rows_removed": 0, "rows_after": n_before}

    mask = df["order_proc_id_coded"].isin(bad_ids)
    df_filtered = df[~mask]
    n_removed = int(mask.sum())
    df_filtered.to_csv(out_path, index=False)

    return {
        "rows_before": n_before,
        "rows_removed": n_removed,
        "rows_after": len(df_filtered),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Remove 'different_patients' conflict rows from all CSVs."
    )
    parser.add_argument("--conflict", default=DEFAULT_CONFLICT,
                        help="Path to the CONFLICT CSV file.")
    parser.add_argument("--input", default=DEFAULT_INPUT,
                        help="Folder containing CSVs to filter.")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help="Folder to write filtered CSVs into.")
    args = parser.parse_args()

    conflict_path = os.path.abspath(args.conflict)
    input_dir     = os.path.abspath(args.input)
    output_dir    = os.path.abspath(args.output)

    # Validate inputs
    if not os.path.isfile(conflict_path):
        print(f"ERROR: Conflict file not found: {conflict_path}")
        sys.exit(1)
    if not os.path.isdir(input_dir):
        print(f"ERROR: Input folder not found: {input_dir}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    bad_ids = load_conflict_ids(conflict_path)

    # Collect CSVs to process — skip CONFLICT_ files and subdirectories
    csv_files = sorted(
        f for f in os.listdir(input_dir)
        if f.lower().endswith(".csv") and not f.startswith("CONFLICT_")
    )

    if not csv_files:
        print(f"No eligible CSV files found in: {input_dir}")
        sys.exit(1)

    print(f"\nProcessing {len(csv_files)} file(s) from: {input_dir}")
    print(f"Writing filtered files to:              {output_dir}\n")

    total_removed = 0

    for fname in csv_files:
        in_path  = os.path.join(input_dir, fname)
        out_path = os.path.join(output_dir, fname)
        summary  = filter_file(in_path, bad_ids, out_path)
        total_removed += summary["rows_removed"]

        status = f"  removed {summary['rows_removed']:>6,} row(s)" if summary["rows_removed"] else "  no rows removed"
        print(f"  {fname:<60s}  {summary['rows_before']:>8,} → {summary['rows_after']:>8,}  |  {status}")

    print(f"\nDone. Total rows removed across all files: {total_removed:,}")


if __name__ == "__main__":
    main()
