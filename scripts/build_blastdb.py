"""Step 3: Download complete Influenza A CDS sequences and build BLAST database.

Downloads ~1.3M sequences from NCBI Entrez in batches of 10,000.
Resulting FASTA is ~2.3 GB; BLAST DB is ~1.5-2 GB.
"""

import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from Bio import Entrez


def run(config, data_dir):
    db_name = config["blast_db_name"]
    entrez_query = config["entrez_query"]
    entrez_email = config["entrez_email"]
    batch_size = config.get("entrez_batch_size", 10000)

    Entrez.email = entrez_email

    db_dir = data_dir / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = db_dir / "genomes.fasta"

    # Check if FASTA already exists and is large enough to skip download
    if fasta_path.exists() and fasta_path.stat().st_size > 1_000_000_000:
        print(f"FASTA already exists ({fasta_path.stat().st_size / 1e9:.1f} GB). Skipping download.")
        print("  Delete it to force re-download.")
        n_seqs = _count_fasta_records(fasta_path)
        print(f"  Sequences in FASTA: {n_seqs:,}")
    else:
        _download_sequences(entrez_query, batch_size, fasta_path)

    # Check if BLAST DB already exists
    db_files = list(db_dir.glob(f"{db_name}.n*"))
    if db_files:
        print(f"BLAST database already exists ({len(db_files)} files). Skipping build.")
        print("  Delete them to force rebuild.")
        return db_dir / db_name

    _build_blast_db(fasta_path, db_dir, db_name)

    return db_dir / db_name


def _download_sequences(entrez_query, batch_size, fasta_path):
    """Download all matching sequences from NCBI Entrez in batches."""
    # First, get total count and WebEnv session
    print(f"Searching NCBI: {entrez_query}")
    handle = Entrez.esearch(
        db="nucleotide",
        term=entrez_query,
        retmax=0,
        usehistory="y",
    )
    record = Entrez.read(handle)
    handle.close()

    total = int(record["Count"])
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    print(f"  Found {total:,} sequences")
    print(f"  Will download in batches of {batch_size:,}")

    n_batches = (total + batch_size - 1) // batch_size
    print(f"  Total batches: {n_batches}")

    # Download in batches, appending to FASTA file
    with open(fasta_path, "w") as fasta_out:
        for i in range(n_batches):
            start = i * batch_size
            attempt = 0
            max_attempts = 5

            while attempt < max_attempts:
                try:
                    print(f"  Batch {i + 1}/{n_batches} (records {start:,}-{min(start + batch_size, total):,})...",
                          end="", flush=True)
                    handle = Entrez.efetch(
                        db="nucleotide",
                        rettype="fasta",
                        retmode="text",
                        retstart=start,
                        retmax=batch_size,
                        webenv=webenv,
                        query_key=query_key,
                    )
                    data = handle.read()
                    handle.close()

                    fasta_out.write(data)
                    print(f" OK ({len(data):,} bytes)")
                    break

                except Exception as e:
                    attempt += 1
                    wait = 10 * attempt
                    print(f" ERROR: {e}")
                    if attempt < max_attempts:
                        print(f"    Retrying in {wait}s (attempt {attempt}/{max_attempts})...")
                        time.sleep(wait)
                    else:
                        print(f"    FAILED after {max_attempts} attempts. Continuing with next batch.")

            # Rate limit: NCBI asks for max 3 requests/second without API key
            time.sleep(0.5)

    file_size = fasta_path.stat().st_size
    n_seqs = _count_fasta_records(fasta_path)
    print(f"\n  FASTA written to {fasta_path}")
    print(f"  File size: {file_size / 1e9:.2f} GB")
    print(f"  Total sequences: {n_seqs:,}")


def _count_fasta_records(fasta_path):
    """Count '>' lines in a FASTA file without loading it into memory."""
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def _build_blast_db(fasta_path, db_dir, db_name):
    """Build BLAST nucleotide database from FASTA file."""
    print("Building BLAST database...")
    print("  (This may take several minutes for a large FASTA)")

    # BLAST+ can't handle paths with spaces, so use a temp dir
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_fasta = os.path.join(tmpdir, "genomes.fasta")

        # Symlink instead of copy for large files
        os.symlink(str(fasta_path.resolve()), tmp_fasta)

        tmp_db = os.path.join(tmpdir, db_name)

        cmd = [
            "makeblastdb",
            "-in", tmp_fasta,
            "-dbtype", "nucl",
            "-out", tmp_db,
            "-parse_seqids",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"STDERR: {result.stderr}")
            raise RuntimeError(f"makeblastdb failed with exit code {result.returncode}")

        print(result.stdout)

        # Copy DB files back to data directory
        for f in os.listdir(tmpdir):
            if f.startswith(db_name) and f != "genomes.fasta":
                shutil.copy2(os.path.join(tmpdir, f), str(db_dir / f))

    print(f"  BLAST database built at {db_dir / db_name}")
