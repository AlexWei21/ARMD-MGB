"""
deidentify_time.py
------------------
Applies per-patient date shifts to every CSV in a target folder.

For each file:
  - Any column whose name matches a known raw date column is shifted and
    renamed <col>_shifted (date only, no time component).
  - Any column already named <col>_shifted (e.g. from a prior pipeline) is
    shifted in place and kept as <col>_shifted.
  - The original raw date columns are dropped after shifting.
  - The date_offset merge column is always dropped from the output.

The shift magnitude comes from date_shift_map.csv, which maps anon_id to a
timedelta string (one consistent random offset per patient, generated in
microbiology_cultures_cohort_query.ipynb).

Usage:
    python deidentify_time.py [--input PATH] [--output PATH] [--shift-map PATH]

Defaults:
    --input      data/new/data/
    --output     data/new/data/deidentified/
    --shift-map  data/new/data/deidentified/date_shift_map.csv
"""

import argparse
import os
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Defaults (relative to the repo root)
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

DEFAULT_INPUT     = os.path.join(REPO_ROOT, "data", "new", "data")
DEFAULT_OUTPUT    = os.path.join(REPO_ROOT, "data", "new", "data", "deidentified")
DEFAULT_SHIFT_MAP = os.path.join(REPO_ROOT, "data", "new", "data", "deidentified", "date_shift_map.csv")

# Known raw date column names (without _shifted suffix).
# Note: 'recorded_resistsant_time' typo preserved to match source data header.
KNOWN_DATE_COLS = {
    "order_time_jittered_utc",
    "medicationOrderDTS",
    "visit_date",
    "recorded_resistsant_time",
    "prior_infecting_organism_recorded_time",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_shift_map(path: str) -> pd.DataFrame:
    """Load date_shift_map.csv and parse date_offset as timedelta."""
    df = pd.read_csv(path)
    if "anon_id" not in df.columns or "date_offset" not in df.columns:
        raise ValueError(f"date_shift_map must have 'anon_id' and 'date_offset' columns: {path}")
    df["date_offset"] = pd.to_timedelta(df["date_offset"])
    return df[["anon_id", "date_offset"]]


def shift_file(in_path: str, out_path: str, shift_map: pd.DataFrame) -> dict:
    """
    Apply date shifts to one CSV file.

    Returns a summary dict with columns shifted and row count.
    """
    df = pd.read_csv(in_path, low_memory=False)

    if "anon_id" not in df.columns:
        # No patient ID to join on — copy unchanged.
        df.to_csv(out_path, index=False)
        return {"raw_shifted": [], "pre_shifted": [], "n_rows": len(df)}

    # Identify which date columns are present
    raw_cols  = [c for c in df.columns if c in KNOWN_DATE_COLS]
    pre_cols  = [c for c in df.columns if c.endswith("_shifted")
                 and c[: -len("_shifted")] in KNOWN_DATE_COLS]

    if not raw_cols and not pre_cols:
        # Nothing to shift — copy unchanged.
        df.to_csv(out_path, index=False)
        return {"raw_shifted": [], "pre_shifted": [], "n_rows": len(df)}

    # Merge offset (left join preserves all rows; unmatched get NaT after shift)
    df = df.merge(shift_map, on="anon_id", how="left")

    # Shift raw date columns → <col>_shifted, drop originals
    for col in raw_cols:
        df[col] = pd.to_datetime(df[col], errors="coerce")
        df[col + "_shifted"] = df[col] + df["date_offset"]
        df.drop(columns=[col], inplace=True)

    # Shift already-shifted columns in place
    for col in pre_cols:
        df[col] = pd.to_datetime(df[col], errors="coerce")
        df[col] = df[col] + df["date_offset"]

    df.drop(columns=["date_offset"], inplace=True)
    df.to_csv(out_path, index=False)

    return {"raw_shifted": raw_cols, "pre_shifted": pre_cols, "n_rows": len(df)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Apply per-patient date shifts to CSV files.")
    parser.add_argument("--input",     default=DEFAULT_INPUT,
                        help="Folder containing CSVs to shift.")
    parser.add_argument("--output",    default=DEFAULT_OUTPUT,
                        help="Folder to write date-shifted CSVs into.")
    parser.add_argument("--shift-map", default=DEFAULT_SHIFT_MAP,
                        help="Path to date_shift_map.csv.")
    args = parser.parse_args()

    input_dir  = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    shift_path = os.path.abspath(args.shift_map)

    # Validate
    if not os.path.isdir(input_dir):
        print(f"ERROR: Input folder not found: {input_dir}")
        sys.exit(1)
    if not os.path.isfile(shift_path):
        print(f"ERROR: Shift map not found: {shift_path}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading shift map from: {shift_path}")
    shift_map = load_shift_map(shift_path)
    print(f"  {len(shift_map):,} patient offsets loaded.\n")

    # Collect CSVs — skip CONFLICT_ files, skip the shift map itself
    shift_map_fname = os.path.basename(shift_path)
    csv_files = sorted(
        f for f in os.listdir(input_dir)
        if f.lower().endswith(".csv")
        and not f.startswith("CONFLICT_")
        and f != shift_map_fname
    )

    if not csv_files:
        print(f"No eligible CSV files found in: {input_dir}")
        sys.exit(1)

    print(f"Processing {len(csv_files)} file(s):")
    print(f"  Input  : {input_dir}")
    print(f"  Output : {output_dir}\n")

    for fname in csv_files:
        in_path  = os.path.join(input_dir, fname)
        out_path = os.path.join(output_dir, fname)
        summary  = shift_file(in_path, out_path, shift_map)

        shifted = summary["raw_shifted"] + summary["pre_shifted"]
        if shifted:
            cols_str = ", ".join(shifted)
            print(f"  [shifted] {fname:<55s}  cols: {cols_str}")
        else:
            print(f"  [copied]  {fname:<55s}  (no date columns)")

    print("\nDone.")


if __name__ == "__main__":
    main()
