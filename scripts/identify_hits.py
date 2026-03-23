"""Step 1: Identify Influenza A hits in precomputed virus data.

Uses species column (s__Alphainfluenzavirus influenzae) as primary filter
to catch all naming conventions (t__Influenza A H3N2, t__H3N2, etc.).
Falls back to subspecies regex patterns if species column is missing.
"""

import re
import pandas as pd


def run(config, results_dir):
    virus_path = config["precomputed_virus_data"]

    print(f"Reading precomputed virus data from {virus_path}")
    df = pd.read_csv(virus_path, sep="\t", compression="gzip", dtype=str)
    print(f"  Total rows: {len(df):,}")

    # Primary: filter by species column
    species_filter = config.get("species_filter")
    if species_filter and "species" in df.columns:
        mask = df["species"] == species_filter
        print(f"  Filtering by species == '{species_filter}'")
    else:
        # Fallback: subspecies regex patterns
        patterns = config["subspecies_patterns"]
        combined = "|".join(f"(?:{p})" for p in patterns)
        mask = df["subspecies"].str.contains(combined, regex=True, na=False)
        print(f"  Filtering by subspecies patterns (species column not available)")

    hits = df.loc[mask].copy()
    print(f"  Influenza A hits: {len(hits):,}")

    if hits.empty:
        print("WARNING: No Influenza A hits found.")

    # Select and output relevant columns
    keep_cols = [
        "sample_ID", "Accession", "Name", "subspecies", "City", "Date",
        "delivery_date", "read_count", "reads_per_million",
    ]
    keep_cols = [c for c in keep_cols if c in hits.columns]
    hits = hits[keep_cols]

    out_path = results_dir / "influenza_hits.tsv"
    hits.to_csv(out_path, sep="\t", index=False)
    print(f"  Written to {out_path}")
    print(f"  Unique samples: {hits['sample_ID'].nunique()}")
    print(f"  Unique subspecies: {hits['subspecies'].unique().tolist()}")

    # Summary by subtype
    if not hits.empty:
        print(f"\n  Hits by EsViritu subtype:")
        for subtype, group in hits.groupby("subspecies"):
            total_reads = group["read_count"].astype(int).sum()
            print(f"    {subtype}: {len(group)} rows, {total_reads:,} reads")

    return hits
