#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import os

# ---------- Helpers ----------
def format_size(bytes_size):
    """Convert bytes to human-readable MB/GB string."""
    gb = bytes_size / (1024 ** 3)
    mb = bytes_size / (1024 ** 2)
    return f"{gb:.2f} GB" if gb >= 1 else f"{mb:.2f} MB"

def is_text(path: Path) -> bool:
    return path.suffix.lower() in {".csv", ".tsv"} or (
        path.suffix.lower() == ".gz" and path.with_suffix("").suffix.lower() in {".csv", ".tsv"}
    )

def is_parquet(path: Path) -> bool:
    return path.suffix.lower() == ".parquet"

def read_table(path: Path) -> pd.DataFrame:
    if is_text(path):
        sep = "," if (path.suffix.lower() == ".csv" or path.with_suffix("").suffix.lower() == ".csv") else "\t"
        return pd.read_csv(path, sep=sep, low_memory=False)
    if is_parquet(path):
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported file type: {path}")

def write_table(df: pd.DataFrame, out_path: Path, like_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if is_text(like_path):
        sep = "," if (like_path.suffix.lower() == ".csv" or like_path.with_suffix("").suffix.lower() == ".csv") else "\t"
        df.to_csv(out_path, index=False, sep=sep)
    elif is_parquet(like_path):
        df.to_parquet(out_path, index=False)
    else:
        raise ValueError(f"Unsupported file type: {like_path}")

def inspect_dataframe(df: pd.DataFrame):
    print(f"  Rows: {len(df):,}")
    print(f"  Columns: {len(df.columns)}")
    print("  Column ranges:")
    for col in df.columns:
        s = df[col]
        dtype = s.dtype
        if pd.api.types.is_numeric_dtype(s):
            print(f"    {col} ({dtype}): min={s.min()}, max={s.max()}")
        elif pd.api.types.is_datetime64_any_dtype(s):
            print(f"    {col} ({dtype}): min={s.min()}, max={s.max()}")
        else:
            # try convert to datetime for inspection
            try:
                s_dt = pd.to_datetime(s, errors='coerce', format='%Y-%m-%d', exact=False)
                if s_dt.notna().any():
                    print(f"    {col} (datetime-like): min={s_dt.min()}, max={s_dt.max()}")
                else:
                    print(f"    {col} ({dtype}): {s.nunique()} unique values")
            except Exception:
                print(f"    {col} ({dtype}): {s.nunique()} unique values")

def deduplicate_rows():
    
    src = Path("deid_time_jittered_data")
    dst = Path("dedup")
    dst.mkdir(exist_ok=True)

    files = [p for p in src.glob("**/*") if p.is_file() and (is_text(p) or is_parquet(p))]
    files.sort()
    
    for f in files:
        rel = f.relative_to(src)
        out_path = dst / rel
        out_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"\nProcessing: {f}")
        df = read_table(f)
        before = len(df)
        df = df.drop_duplicates()
        after = len(df)
        print(f"  Deduplicated: {before:,} → {after:,} (removed {before - after:,})")
        inspect_dataframe(df)

        write_table(df, out_path, f)
        print(f"  Saved to {out_path}")

        before_size = f.stat().st_size
        after_size = out_path.stat().st_size
        print(f"  File size before: {format_size(before_size)} -> after: {format_size(after_size)}")

    print("\nAll files processed successfully.")

def main():

    deduplicate_rows()
    
if __name__ == "__main__":
    main()
