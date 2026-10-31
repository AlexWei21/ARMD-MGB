"""
filter_multi_encounter_cultures.py
----------------------------------
Removes 'same_patient_diff_encounter' conflicts from an already-deidentified
data pull: cultures (order_proc_id_coded) that appear with more than one
pat_enc_csn_id_coded in the cohort file.

These arise in microbiology_cultures_cohort_query.ipynb when a culture's
collection time falls inside two or more overlapping encounter windows
(e.g. ED + Inpatient), so every row of the culture is duplicated once per
encounter. Those orders are removed from every CSV in the input folder that
has an order_proc_id_coded column.

Run after filter_deid_conflicts.py (which removes 'different_patients'
conflicts), so each order_proc_id_coded belongs to a single patient.

Usage:
    python filter_multi_encounter_cultures.py --input <folder> --output <folder>
        [--cohort microbiology_cohort.csv]
"""

import argparse
import os
import sys

import pandas as pd

DEFAULT_COHORT_NAME = "microbiology_cohort.csv"
CHUNKSIZE = 2_000_000


def find_multi_encounter_orders(cohort_path: str) -> set:
    """Return order_proc_id_coded values that map to >1 pat_enc_csn_id_coded."""
    df = pd.read_csv(cohort_path, usecols=["order_proc_id_coded", "pat_enc_csn_id_coded"])
    counts = df.drop_duplicates().groupby("order_proc_id_coded")["pat_enc_csn_id_coded"].nunique()
    bad = set(counts[counts > 1].index)
    print(f"Detected {len(bad):,} order_proc_id_coded values linked to >1 encounter.")
    return bad


def filter_file(path: str, bad_ids: set, out_path: str) -> dict:
    # Stream in chunks: some files (e.g. comorbidity) are several GB.
    n_before = n_removed = 0
    header = True
    for chunk in pd.read_csv(path, low_memory=False, chunksize=CHUNKSIZE):
        n_before += len(chunk)
        if "order_proc_id_coded" in chunk.columns:
            mask = chunk["order_proc_id_coded"].isin(bad_ids)
            n_removed += int(mask.sum())
            chunk = chunk[~mask]
        chunk.to_csv(out_path, index=False, mode="w" if header else "a", header=header)
        header = False
    return {
        "rows_before": n_before,
        "rows_removed": n_removed,
        "rows_after": n_before - n_removed,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input folder with deidentified CSVs.")
    parser.add_argument("--output", required=True, help="Output folder for filtered CSVs.")
    parser.add_argument("--cohort", default=DEFAULT_COHORT_NAME,
                        help="Cohort file name (used to detect multi-encounter orders).")
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
    if input_dir == output_dir:
        print("ERROR: Output folder must differ from input folder.")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    bad_ids = find_multi_encounter_orders(cohort_path)

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
        print(f"  {fname:<65s}  {summary['rows_before']:>10,} → {summary['rows_after']:>10,}  |  {tag}", flush=True)

    print(f"\nDone. Total rows removed across all files: {total_removed:,}")


if __name__ == "__main__":
    main()
