"""Step 2: Download BAMs from GCS and extract Influenza A-aligned reads."""

import os
import tempfile
import pandas as pd
import pysam
from google.cloud import storage


def _get_gcs_client(config):
    sa_path = config["gcs_service_account"]
    return storage.Client.from_service_account_json(sa_path)


def _delivery_date_to_folder(delivery_date):
    """Convert delivery_date (YYYY-MM-DD) to YYYYMMDD folder prefix."""
    return delivery_date.replace("-", "")


def _find_bam_path(client, bucket_name, prefix, delivery_date, sample_id):
    """Find the BAM file path in GCS for a given sample and delivery."""
    folder = _delivery_date_to_folder(delivery_date)
    bam_name = f"{sample_id}.third.filt.sorted.bam"
    bam_path = f"{prefix}/{folder}_esviritu_outputs/{bam_name}"

    bucket = client.bucket(bucket_name)
    blob = bucket.blob(bam_path)
    if blob.exists():
        return bam_path
    return None


def _download_blob(client, bucket_name, blob_path, local_path):
    """Download a blob from GCS to local path."""
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_path)
    blob.download_to_filename(local_path)


def run(config, results_dir):
    hits_path = results_dir / "influenza_hits.tsv"
    hits = pd.read_csv(hits_path, sep="\t", dtype=str)

    if hits.empty:
        print("No Influenza A hits to extract reads from.")
        (results_dir / "influenza_reads.fasta").write_text("")
        pd.DataFrame().to_csv(results_dir / "read_metadata.tsv", sep="\t", index=False)
        return

    client = _get_gcs_client(config)
    bucket_name = config["gcs_bucket"]
    prefix = config["gcs_prefix"]

    # Group by (sample_ID, delivery_date) to avoid duplicate BAM downloads
    grouped = hits.groupby(["sample_ID", "delivery_date"])

    # Collect accessions per (sample, delivery) for targeted read extraction
    sample_accessions = {}
    sample_meta = {}
    for (sample_id, delivery_date), group in grouped:
        key = (sample_id, delivery_date)
        sample_accessions[key] = set(group["Accession"].dropna().unique())
        row = group.iloc[0]
        sample_meta[key] = {
            "City": row.get("City", ""),
            "Date": row.get("Date", ""),
            "subspecies_list": ", ".join(group["subspecies"].unique()),
        }

    fasta_records = []
    metadata_rows = []
    total_reads = 0

    for (sample_id, delivery_date), accessions in sample_accessions.items():
        meta = sample_meta[(sample_id, delivery_date)]
        print(f"Processing {sample_id} (delivery {delivery_date})...")

        bam_path = _find_bam_path(client, bucket_name, prefix, delivery_date, sample_id)
        if bam_path is None:
            print(f"  WARNING: BAM not found for {sample_id} in delivery {delivery_date}")
            continue

        bai_path = bam_path + ".bai"

        with tempfile.TemporaryDirectory() as tmpdir:
            local_bam = os.path.join(tmpdir, f"{sample_id}.bam")
            local_bai = os.path.join(tmpdir, f"{sample_id}.bam.bai")

            try:
                print(f"  Downloading BAM...")
                _download_blob(client, bucket_name, bam_path, local_bam)
            except Exception as e:
                print(f"  ERROR downloading BAM: {e}")
                continue

            try:
                _download_blob(client, bucket_name, bai_path, local_bai)
            except Exception as e:
                print(f"  WARNING: BAI not found, indexing locally...")
                try:
                    pysam.index(local_bam)
                except Exception as e2:
                    print(f"  ERROR indexing BAM: {e2}")
                    continue

            try:
                bam = pysam.AlignmentFile(local_bam, "rb")
            except Exception as e:
                print(f"  ERROR opening BAM: {e}")
                continue

            # Get reference names from BAM header
            ref_names = set(bam.references)
            target_refs = accessions & ref_names

            if not target_refs:
                print(f"  WARNING: Accessions {accessions} not in BAM references")
                # Try partial matching (strip version number)
                for acc in accessions:
                    acc_base = acc.split(".")[0]
                    for ref in ref_names:
                        if acc_base in ref:
                            target_refs.add(ref)
                if target_refs:
                    print(f"  Found partial matches: {target_refs}")
                else:
                    bam.close()
                    continue

            sample_read_count = 0
            for ref in target_refs:
                for read in bam.fetch(ref):
                    if read.is_unmapped or read.query_sequence is None:
                        continue

                    read_name = read.query_name
                    seq = read.query_sequence
                    city = meta["City"]
                    date = meta["Date"]
                    subspecies = meta["subspecies_list"]

                    uid = f"read_{total_reads + sample_read_count}"
                    header = f">{uid} {read_name}|{sample_id}|{city}|{date}|{subspecies}|{ref}"
                    fasta_records.append(f"{header}\n{seq}")

                    metadata_rows.append({
                        "read_uid": uid,
                        "read_id": read_name,
                        "sample_ID": sample_id,
                        "City": city,
                        "Date": date,
                        "delivery_date": delivery_date,
                        "subspecies": subspecies,
                        "accession": ref,
                        "query_length": len(seq),
                        "read_sequence": seq,
                    })
                    sample_read_count += 1

            bam.close()
            total_reads += sample_read_count
            print(f"  Extracted {sample_read_count} reads from {len(target_refs)} reference(s)")

    # Write FASTA
    fasta_path = results_dir / "influenza_reads.fasta"
    fasta_path.write_text("\n".join(fasta_records) + "\n" if fasta_records else "")
    print(f"\nTotal reads extracted: {total_reads}")
    print(f"FASTA written to {fasta_path}")

    # Write metadata TSV
    meta_df = pd.DataFrame(metadata_rows)
    meta_path = results_dir / "read_metadata.tsv"
    meta_df.to_csv(meta_path, sep="\t", index=False)
    print(f"Read metadata written to {meta_path}")

    return meta_df
