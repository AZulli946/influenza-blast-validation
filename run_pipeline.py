#!/usr/bin/env python3
"""Influenza A BLAST Validation Pipeline — main orchestrator.

Usage:
    python run_pipeline.py [--config config.yaml] [--step STEP]

Steps: identify, extract, blastdb, blast, subtype, all (default)
"""

import argparse
import sys
import time
from pathlib import Path

import yaml


def load_config(config_path):
    """Load config, resolving relative paths against the project directory."""
    project_dir = Path(__file__).resolve().parent
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Resolve relative paths
    for key in ("precomputed_virus_data", "gcs_service_account"):
        if key in config and not Path(config[key]).is_absolute():
            config[key] = str(project_dir / config[key])

    return config, project_dir


def main():
    parser = argparse.ArgumentParser(description="Influenza A BLAST Validation Pipeline")
    parser.add_argument(
        "--config", default="config.yaml",
        help="Path to config.yaml (default: config.yaml in script directory)",
    )
    parser.add_argument(
        "--step", default="all",
        choices=["identify", "extract", "blastdb", "blast", "subtype", "all"],
        help="Run a specific step (default: all)",
    )
    args = parser.parse_args()

    # Resolve config path relative to script directory
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = Path(__file__).resolve().parent / config_path

    config, project_dir = load_config(config_path)

    results_dir = project_dir / "results"
    data_dir = project_dir / "data"
    results_dir.mkdir(exist_ok=True)
    data_dir.mkdir(exist_ok=True)

    step = args.step
    t0 = time.time()

    # Step 1: Identify influenza A hits
    if step in ("identify", "all"):
        print("=" * 60)
        print("STEP 1: Identify Influenza A hits")
        print("=" * 60)
        from scripts.identify_hits import run as identify
        identify(config, results_dir)
        print()

    # Step 2: Extract reads from BAMs
    if step in ("extract", "all"):
        print("=" * 60)
        print("STEP 2: Extract reads from GCS BAMs")
        print("=" * 60)
        from scripts.extract_reads import run as extract
        extract(config, results_dir)
        print()

    # Step 3: Build BLAST database
    if step in ("blastdb", "all"):
        print("=" * 60)
        print("STEP 3: Build BLAST database")
        print("=" * 60)
        from scripts.build_blastdb import run as build_db
        build_db(config, data_dir)
        print()

    # Step 4: Run BLAST validation
    if step in ("blast", "all"):
        print("=" * 60)
        print("STEP 4: BLAST validation")
        print("=" * 60)
        from scripts.blast_validate import run as blast
        blast(config, data_dir, results_dir)
        print()

    # Step 5: Assign subtypes
    if step in ("subtype", "all"):
        print("=" * 60)
        print("STEP 5: Assign subtypes")
        print("=" * 60)
        from scripts.assign_subtypes import run as subtype
        subtype(config, project_dir, results_dir)
        print()

    elapsed = time.time() - t0
    print(f"Pipeline complete in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
