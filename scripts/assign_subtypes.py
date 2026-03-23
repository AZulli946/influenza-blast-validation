"""Step 5: Assign segment identity and subtype from BLAST hits.

For each read:
1. Parse segment (HA, NA, PB2, PB1, PA, NP, M, NS) from BLAST hit titles
2. Parse subtype (H1N1, H3N2, etc.) from BLAST hit titles
3. Majority-vote subtype using preferred hits (>=98% identity, non-terminal) first;
   falls back to all hits if no preferred hits exist for a read.
4. Confidence = count(majority subtype) / count(total hits with parseable subtype)
5. identity_tier = "high" if assignment used preferred hits, "low" if fallback

Output: one row per read with subtype assignment, confidence, segment, and tier.
All individual BLAST hits are also preserved in a separate file.
"""

import re
import pandas as pd
from collections import Counter


# Regex for parsing strain names like A/swine/Minnesota/A029807704/2026(H1N2)
_SUBTYPE_RE = re.compile(r'\(H(\d+)N(\d+)\)')
_H_ONLY_RE = re.compile(r'\(H(\d+)\)')
_MIXED_RE = re.compile(r'\bmixed\b', re.IGNORECASE)

# Segment identification from BLAST hit titles
_SEGMENT_NUMBER_RE = re.compile(r'segment\s+(\d)')
_SEGMENT_NAME_RE = re.compile(
    r'\b(hemagglutinin|neuraminidase|polymerase\s+(?:basic\s+(?:protein\s+)?[12]|acidic)|'
    r'nucleoprotein|matrix|nonstructural|nuclear\s+export)\b',
    re.IGNORECASE,
)
_SEGMENT_ABBREV_RE = re.compile(r'\b(HA|NA|PB2|PB1|PA|NP|NS|NEP)\b')

SEGMENT_NUM_TO_NAME = {
    "1": "PB2",
    "2": "PB1",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "M",
    "8": "NS",
}

SEGMENT_FULLNAME_TO_ABBREV = {
    "hemagglutinin": "HA",
    "neuraminidase": "NA",
    "nucleoprotein": "NP",
    "matrix": "M",
    "nonstructural": "NS",
    "nuclear export": "NS",
}


def _parse_subtype(stitle):
    """Extract subtype (e.g., 'H3N2') from BLAST subject title."""
    if pd.isna(stitle):
        return None
    if _MIXED_RE.search(stitle):
        return None

    m = _SUBTYPE_RE.search(stitle)
    if m:
        return f"H{m.group(1)}N{m.group(2)}"

    m = _H_ONLY_RE.search(stitle)
    if m:
        return f"H{m.group(1)}"

    return None


def _parse_segment(stitle):
    """Extract segment identity from BLAST subject title."""
    if pd.isna(stitle):
        return None

    m = _SEGMENT_NUMBER_RE.search(stitle)
    if m:
        return SEGMENT_NUM_TO_NAME.get(m.group(1))

    m = _SEGMENT_ABBREV_RE.search(stitle)
    if m:
        abbrev = m.group(1)
        if abbrev == "NEP":
            return "NS"
        return abbrev

    m = _SEGMENT_NAME_RE.search(stitle)
    if m:
        name = m.group(1).lower()
        if "polymerase" in name:
            if "basic" in name:
                if "1" in name:
                    return "PB1"
                return "PB2"
            return "PA"
        return SEGMENT_FULLNAME_TO_ABBREV.get(name)

    return None


def _majority_vote(values):
    """Return (majority_value, confidence) from a list of values."""
    if not values:
        return None, 0.0

    counts = Counter(values)
    majority_val, majority_count = counts.most_common(1)[0]
    confidence = majority_count / len(values)
    return majority_val, confidence


def run(config, project_dir, results_dir):
    merged_path = results_dir / "blast_merged.tsv"

    if not merged_path.exists():
        print("No blast_merged.tsv found. Skipping subtype assignment.")
        return

    min_pident = config.get("min_percent_identity", 98.0)

    df = pd.read_csv(merged_path, sep="\t", dtype=str)

    if df.empty:
        print("No data to process.")
        df.to_csv(results_dir / "blast_validated.tsv", sep="\t", index=False)
        return

    # Parse segment and subtype for every BLAST hit row
    print("Parsing segment and subtype from BLAST hit titles...")
    df["blast_segment"] = df["stitle"].apply(_parse_segment)
    df["blast_subtype"] = df["stitle"].apply(_parse_subtype)

    # Ensure flag columns are boolean (they come in as strings from TSV)
    for col in ("above_identity_threshold", "in_terminal_region", "preferred_hit"):
        if col in df.columns:
            df[col] = df[col].map({"True": True, "False": False, True: True, False: False})

    # Report parsing success rates
    has_hit = df["stitle"].notna()
    n_hits = has_hit.sum()
    n_segment_parsed = df.loc[has_hit, "blast_segment"].notna().sum()
    n_subtype_parsed = df.loc[has_hit, "blast_subtype"].notna().sum()
    print(f"  BLAST hit rows: {n_hits:,}")
    print(f"  Segment parsed: {n_segment_parsed:,} ({100 * n_segment_parsed / max(n_hits, 1):.1f}%)")
    print(f"  Subtype parsed: {n_subtype_parsed:,} ({100 * n_subtype_parsed / max(n_hits, 1):.1f}%)")

    # Save all hits with parsed fields (many rows per read)
    all_hits_path = results_dir / "blast_all_hits.tsv"
    df.to_csv(all_hits_path, sep="\t", index=False)
    print(f"  All hits with parsed fields written to {all_hits_path}")

    # Aggregate per read: majority-vote segment and subtype
    print("\nAssigning per-read segment and subtype by majority vote...")
    print(f"  Using preferred hits (>={min_pident}% identity, non-terminal) when available")
    print(f"  Falling back to all hits when no preferred hits exist")

    meta_cols = ["read_uid", "read_id", "sample_ID", "City", "Date",
                 "delivery_date", "subspecies", "accession", "query_length", "read_sequence"]
    meta_cols = [c for c in meta_cols if c in df.columns]

    read_groups = df.groupby("read_uid")

    read_results = []
    n_high_tier = 0
    n_low_tier = 0
    n_no_hits = 0

    for read_uid, group in read_groups:
        meta = group.iloc[0][meta_cols].to_dict()

        hit_rows = group[group["stitle"].notna()]
        n_hits_read = len(hit_rows)

        if n_hits_read == 0:
            meta.update({
                "n_blast_hits": 0,
                "n_preferred_hits": 0,
                "identity_tier": "none",
                "assigned_segment": None,
                "segment_confidence": 0.0,
                "assigned_subtype": None,
                "subtype_confidence": 0.0,
                "top_blast_hit": None,
                "top_pident": None,
                "top_evalue": None,
                "top_bitscore": None,
            })
            read_results.append(meta)
            n_no_hits += 1
            continue

        # Determine which hits to use for assignment
        preferred = hit_rows[hit_rows["preferred_hit"] == True] if "preferred_hit" in hit_rows.columns else pd.DataFrame()
        n_preferred = len(preferred)

        if n_preferred > 0:
            assignment_rows = preferred
            tier = "high"
            n_high_tier += 1
        else:
            assignment_rows = hit_rows
            tier = "low"
            n_low_tier += 1

        # Majority vote on segment (from assignment rows)
        segments = assignment_rows["blast_segment"].dropna().tolist()
        assigned_segment, segment_conf = _majority_vote(segments)

        # Majority vote on subtype (from assignment rows)
        subtypes = assignment_rows["blast_subtype"].dropna().tolist()
        assigned_subtype, subtype_conf = _majority_vote(subtypes)

        # Top hit by bitscore (from ALL hits, not just preferred)
        hit_rows_numeric = hit_rows.copy()
        hit_rows_numeric["bitscore"] = pd.to_numeric(hit_rows_numeric["bitscore"], errors="coerce")
        top_hit = hit_rows_numeric.loc[hit_rows_numeric["bitscore"].idxmax()]

        meta.update({
            "n_blast_hits": n_hits_read,
            "n_preferred_hits": n_preferred,
            "identity_tier": tier,
            "assigned_segment": assigned_segment,
            "segment_confidence": round(segment_conf, 4),
            "assigned_subtype": assigned_subtype,
            "subtype_confidence": round(subtype_conf, 4),
            "top_blast_hit": top_hit.get("stitle"),
            "top_pident": top_hit.get("pident"),
            "top_evalue": top_hit.get("evalue"),
            "top_bitscore": top_hit.get("bitscore"),
        })
        read_results.append(meta)

    result_df = pd.DataFrame(read_results)

    # Reorder columns
    leading_cols = [
        "read_id", "sample_ID", "City", "Date", "delivery_date",
        "subspecies", "accession",
        "assigned_segment", "segment_confidence",
        "assigned_subtype", "subtype_confidence",
        "identity_tier", "n_blast_hits", "n_preferred_hits",
        "top_blast_hit", "top_pident", "top_evalue", "top_bitscore",
        "query_length", "read_sequence",
    ]
    leading_cols = [c for c in leading_cols if c in result_df.columns]
    extra_cols = [c for c in result_df.columns if c not in leading_cols]
    result_df = result_df[leading_cols + extra_cols]

    out_path = results_dir / "blast_validated.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)

    # Summary statistics
    n_total = len(result_df)
    print(f"\nFinal results: {n_total:,} reads")
    print(f"  High-confidence tier (>={min_pident}% identity): {n_high_tier:,}")
    print(f"  Low-confidence tier (<{min_pident}% identity only): {n_low_tier:,}")
    print(f"  No BLAST hits: {n_no_hits:,}")

    n_with_hits = n_high_tier + n_low_tier
    if n_with_hits > 0:
        validated = result_df[result_df["n_blast_hits"].astype(int) > 0]

        print(f"\n  Segment distribution:")
        seg_counts = validated["assigned_segment"].value_counts(dropna=False)
        for seg, count in seg_counts.items():
            seg_label = seg if pd.notna(seg) else "Unknown"
            print(f"    {seg_label}: {count:,}")

        print(f"\n  Subtype distribution:")
        sub_counts = validated["assigned_subtype"].value_counts(dropna=False)
        for sub, count in sub_counts.items():
            sub_label = sub if pd.notna(sub) else "Unknown"
            print(f"    {sub_label}: {count:,}")

        print(f"\n  Mean segment confidence: {validated['segment_confidence'].astype(float).mean():.3f}")
        print(f"  Mean subtype confidence: {validated['subtype_confidence'].astype(float).mean():.3f}")

        # Tier breakdown by subtype
        print(f"\n  Subtype by identity tier:")
        tier_sub = pd.crosstab(
            validated["assigned_subtype"].fillna("Unknown"),
            validated["identity_tier"],
            margins=True,
        )
        print(tier_sub.to_string())

        # Cross-tab: segment x subtype
        print(f"\n  Segment x Subtype cross-tabulation:")
        xtab = pd.crosstab(
            validated["assigned_segment"].fillna("Unknown"),
            validated["assigned_subtype"].fillna("Unknown"),
            margins=True,
        )
        print(xtab.to_string())

    print(f"\nWritten to {out_path}")

    return result_df
