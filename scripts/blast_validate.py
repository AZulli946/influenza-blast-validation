"""Step 4: Run BLAST validation with top 100 hits, flag by identity and terminal position.

Key differences from poliovirus pipeline:
- max_target_seqs 100 (not 1) — multiple hits per read
- All hits retained; flagged with above_identity_threshold and in_terminal_region
- 98% identity threshold used downstream for subtype assignment (preferred tier)
- Merge produces many-to-one (many hits per read)
"""

import os
import shutil
import subprocess
import tempfile
import pandas as pd


BLAST_COLUMNS = [
    "qseqid", "sseqid", "stitle", "pident", "length",
    "qlen", "sstart", "send", "slen", "evalue", "bitscore",
]


def run(config, data_dir, results_dir):
    db_name = config["blast_db_name"]
    query_path = results_dir / "influenza_reads.fasta"
    blast_out = results_dir / "blast_results_raw.tsv"
    meta_path = results_dir / "read_metadata.tsv"

    if not query_path.exists() or query_path.stat().st_size == 0:
        print("No reads to BLAST (empty FASTA). Skipping.")
        pd.DataFrame().to_csv(results_dir / "blast_validated.tsv", sep="\t", index=False)
        return

    db_files = list(data_dir.glob(f"blastdb/{db_name}.*"))
    if not db_files:
        raise FileNotFoundError("BLAST database not found. Run build_blastdb first.")

    max_target_seqs = config.get("max_target_seqs", 100)
    evalue = config.get("evalue", 1e-5)
    num_threads = config.get("num_threads", 4)
    min_pident = config.get("min_percent_identity", 98.0)
    terminal_bp = config.get("terminal_filter_bp", 150)

    outfmt = "6 " + " ".join(BLAST_COLUMNS)

    # Only re-run BLAST if raw results don't exist yet
    if blast_out.exists() and blast_out.stat().st_size > 0:
        print(f"Raw BLAST results already exist ({blast_out.stat().st_size / 1e6:.0f} MB). Skipping BLAST.")
        print("  Delete blast_results_raw.tsv to force re-run.")
    else:
        # BLAST+ can't handle paths with spaces
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_query = os.path.join(tmpdir, "query.fasta")
            tmp_out = os.path.join(tmpdir, "blast_results.tsv")

            shutil.copy2(str(query_path), tmp_query)

            # Symlink BLAST DB files to tmpdir (avoid copying multi-GB DB)
            db_dir = data_dir / "blastdb"
            for f in db_dir.iterdir():
                if f.name.startswith(db_name):
                    os.symlink(str(f.resolve()), os.path.join(tmpdir, f.name))

            tmp_db = os.path.join(tmpdir, db_name)

            cmd = [
                "blastn",
                "-query", tmp_query,
                "-db", tmp_db,
                "-outfmt", outfmt,
                "-max_target_seqs", str(max_target_seqs),
                "-evalue", str(evalue),
                "-num_threads", str(num_threads),
                "-out", tmp_out,
            ]

            print(f"Running BLAST ({max_target_seqs} max targets, evalue {evalue}, {num_threads} threads)...")
            print(f"  This may take 1-4 hours for ~57k reads against ~1.3M subjects.")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"STDERR: {result.stderr}")
                raise RuntimeError(f"blastn failed with exit code {result.returncode}")

            shutil.copy2(tmp_out, str(blast_out))

    # Parse raw BLAST results
    if blast_out.stat().st_size == 0:
        print("WARNING: BLAST produced no hits.")
        blast_df = pd.DataFrame(columns=BLAST_COLUMNS)
    else:
        blast_df = pd.read_csv(blast_out, sep="\t", header=None, names=BLAST_COLUMNS)

    n_raw = len(blast_df)
    n_queries_raw = blast_df["qseqid"].nunique()
    print(f"  Raw BLAST hits: {n_raw:,} across {n_queries_raw:,} reads")

    # Convert numeric columns
    blast_df["pident"] = pd.to_numeric(blast_df["pident"], errors="coerce")
    blast_df["sstart"] = pd.to_numeric(blast_df["sstart"], errors="coerce")
    blast_df["send"] = pd.to_numeric(blast_df["send"], errors="coerce")
    blast_df["slen"] = pd.to_numeric(blast_df["slen"], errors="coerce")

    # Flag 1: identity threshold (do NOT remove — just flag)
    blast_df["above_identity_threshold"] = blast_df["pident"] >= min_pident
    n_above = blast_df["above_identity_threshold"].sum()
    n_below = len(blast_df) - n_above
    print(f"  Hits >= {min_pident}% identity: {n_above:,}")
    print(f"  Hits <  {min_pident}% identity: {n_below:,}")

    # Flag 2: terminal region (do NOT remove — just flag)
    def _in_terminal_region(row):
        """Check if alignment is entirely within terminal regions of subject."""
        sstart, send, slen = row["sstart"], row["send"], row["slen"]
        if pd.isna(sstart) or pd.isna(send) or pd.isna(slen):
            return False

        aln_start = min(sstart, send)
        aln_end = max(sstart, send)

        if aln_end <= terminal_bp:
            return True
        if aln_start >= (slen - terminal_bp + 1):
            return True

        return False

    blast_df["in_terminal_region"] = blast_df.apply(_in_terminal_region, axis=1)
    n_terminal = blast_df["in_terminal_region"].sum()
    print(f"  Hits in terminal {terminal_bp}bp: {n_terminal:,}")

    # Preferred hits: above threshold AND not in terminal region
    blast_df["preferred_hit"] = blast_df["above_identity_threshold"] & ~blast_df["in_terminal_region"]
    n_preferred = blast_df["preferred_hit"].sum()
    n_preferred_reads = blast_df.loc[blast_df["preferred_hit"], "qseqid"].nunique()
    print(f"  Preferred hits (>={min_pident}% identity, non-terminal): {n_preferred:,} across {n_preferred_reads:,} reads")

    # Reads with BLAST hits but none preferred
    reads_with_hits = set(blast_df["qseqid"].unique())
    reads_with_preferred = set(blast_df.loc[blast_df["preferred_hit"], "qseqid"].unique())
    reads_low_only = reads_with_hits - reads_with_preferred
    print(f"  Reads with only sub-threshold / terminal hits: {len(reads_low_only):,}")

    # Save all flagged BLAST results
    blast_flagged_path = results_dir / "blast_results_flagged.tsv"
    blast_df.to_csv(blast_flagged_path, sep="\t", index=False)
    print(f"  Flagged BLAST results written to {blast_flagged_path}")

    # Load read metadata and merge
    meta_df = pd.read_csv(meta_path, sep="\t", dtype=str)

    # Left join: keep ALL reads, including those with no BLAST hits at all
    merged = meta_df.merge(blast_df, left_on="read_uid", right_on="qseqid", how="left")
    if "qseqid" in merged.columns:
        merged.drop(columns=["qseqid"], inplace=True)

    n_no_hit = merged.loc[merged["sseqid"].isna(), "read_uid"].nunique()
    n_with_hit = merged.loc[merged["sseqid"].notna(), "read_uid"].nunique()
    print(f"\n  Reads with any BLAST hit: {n_with_hit:,}")
    print(f"  Reads without any BLAST hit: {n_no_hit:,}")
    print(f"  Total rows in merged output: {len(merged):,} (multiple hits per read)")

    merged.to_csv(results_dir / "blast_merged.tsv", sep="\t", index=False)
    print(f"  Merged results written to {results_dir / 'blast_merged.tsv'}")

    return merged
