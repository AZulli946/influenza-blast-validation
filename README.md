# Influenza A BLAST Validation Pipeline

Extracts Influenza A virus-classified reads from EsViritu wastewater metagenomics BAM files, validates them via BLAST against all complete IAV coding sequences in NCBI (~1.33M sequences), assigns segment identity and subtype by majority vote across top BLAST hits, and reports per-read confidence scores.

GitHub: https://github.com/AZulli946/influenza-blast-validation

## Quick Start

```bash
cd influenza-blast-validation
pip install -r requirements.txt    # pandas, pysam, google-cloud-storage, biopython, pyyaml
python run_pipeline.py             # runs all 5 steps (uses config.yaml)
```

Run individual steps or with alternate config:

```bash
python run_pipeline.py --step identify   # step 1 only
python run_pipeline.py --step blastdb    # step 3 only (can run in parallel with step 2)
python run_pipeline.py --config config_first_pass.yaml  # use first-pass BAM config
```

## Pipeline Steps

| Step | Script | Description | Runtime |
|------|--------|-------------|---------|
| 1 | `identify_hits.py` | Filter precomputed dashboard data for Influenza A hits | ~2s |
| 2 | `extract_reads.py` | Download BAMs from GCS, extract IAV-aligned reads | ~10 min |
| 3 | `build_blastdb.py` | Download ~1.33M IAV complete CDS from NCBI Entrez, build BLAST DB | ~2 hrs (one-time) |
| 4 | `blast_validate.py` | Run blastn (100 max targets), flag by identity threshold and terminal position | ~75 min |
| 5 | `assign_subtypes.py` | Parse segment/subtype from hits, majority-vote assignment with confidence | ~1.5 min |

Steps 2 and 3 are independent and can be run in parallel. Steps 4 and 5 depend on prior steps completing.

### Step 1: Identify Hits

Reads `precomputed_virus.tsv.gz` from the metagenomics dashboard and filters for Influenza A hits.

**Primary filter**: Species column (`s__Alphainfluenzavirus influenzae`). This catches all EsViritu naming conventions — both `t__Influenza A H3N2` and bare `t__H3N2` labels — because the species column is consistent regardless of subspecies naming.

**Fallback filter**: Subspecies regex patterns (`^t__Influenza A`, `^t__H\d+N\d+`) if the species column is missing. The fallback uses both patterns to capture both naming conventions. Importantly, bare `t__H1`, `t__H2`, `t__H3` (without N) are **Rotavirus**, not Influenza — the species column filter avoids this false-match problem.

**Third-pass results**: 1,954 hits across 310 unique samples.

**Output**: `results/influenza_hits.tsv`

### Step 2: Extract Reads

For each hit, downloads the corresponding BAM file from GCS (`gs://wwscanseq/viralmetagenomics/reports/`), extracts aligned reads by accession number, and writes sequences + metadata.

The BAM suffix is configurable via `bam_suffix` in the config file. The default is `.third.filt.sorted.bam` (third-pass alignment). Set to `.first.filt.sorted.bam` to use first-pass BAMs (see [Multi-Pass BAM Support](#multi-pass-bam-support)).

**Third-pass results**: 57,914 reads extracted from 310 samples.

**Outputs**:
- `results/influenza_reads.fasta` — query sequences for BLAST
- `results/read_metadata.tsv` — per-read metadata (sample, site, date, EsViritu classification)

### Step 3: Build BLAST Database

Downloads all Influenza A virus complete CDS sequences from NCBI Entrez in batches of 10,000 using the search history (WebEnv) API. Builds a local BLAST nucleotide database.

- **Query**: `"Influenza A virus"[Organism] AND "complete cds"[Title]`
- **Database size**: ~1.33M sequences, ~2.4 GB FASTA, ~1.5-2 GB BLAST DB
- Skips download if FASTA already exists (>1 GB); skips DB build if DB files already exist
- Delete `data/blastdb/genomes.fasta` or DB files to force rebuild
- Retry logic: up to 5 attempts per batch with exponential backoff
- Uses symlinks instead of file copies for large FASTA files

**Why Entrez over NCBI IVR FTP**: The Influenza Virus Resource FTP database is frozen at October 2020, missing 6 years of sequences including critical H5N1 clade 2.3.4.4b. Entrez provides the full, current NCBI collection.

**Output**: `data/blastdb/influenza_a_complete.*`

### Step 4: BLAST Validation

Runs `blastn` with up to 100 target sequences per read. **All hits are retained and flagged** — nothing is filtered out:

| Flag | Column | Meaning |
|------|--------|---------|
| Identity threshold | `above_identity_threshold` | True if percent identity >= 98% |
| Terminal region | `in_terminal_region` | True if alignment falls entirely within the first or last 150 bp of the subject sequence |
| Preferred hit | `preferred_hit` | True if above threshold AND not in terminal region |

The 150 bp terminal filter addresses conserved UTR/signal peptide/transmembrane regions at segment ends that can produce cross-subtype BLAST matches. IAV UTRs are 20-46 bp, with an additional ~100 bp of moderately conserved coding sequence beyond them. The 150 bp threshold is a compromise: aggressive enough to catch the conserved region but not so aggressive as to remove too many true hits on shorter segments like M (1,027 bp) and NS (890 bp).

Skips the BLAST search if `blast_results_raw.tsv` already exists. Delete it to force re-run. Uses symlinks for BLAST DB files to avoid copying multi-GB files.

**Outputs**:
- `results/blast_results_raw.tsv` — raw BLAST output (tab-separated, no header)
- `results/blast_results_flagged.tsv` — all hits with flag columns
- `results/blast_merged.tsv` — hits left-joined to read metadata (many rows per read)

### Step 5: Assign Subtypes

Parses segment and subtype from each BLAST hit's subject title, then assigns per-read values by majority vote.

**Segment identification** — parsed from BLAST subject titles using (in priority order):
1. Segment number: `"segment 4"` -> HA
2. Gene abbreviation: `"(HA)"`, `"(NA)"`, `"(PB2)"`
3. Full gene name: `"hemagglutinin"`, `"neuraminidase"`

**Subtype identification** — parsed from strain name parenthetical:
- `A/swine/Minnesota/.../2026(H1N2)` -> H1N2
- `A/duck/.../2024(H5)` -> H5 (NA not annotated)
- Mixed subtypes are excluded

**Tiered assignment**:
- **High tier**: Read has preferred hits (>=98% identity, non-terminal). Majority vote uses only preferred hits.
- **Low tier**: Read has BLAST hits but none are preferred. Majority vote uses all hits. These reads are still assigned a subtype but should be interpreted with caution.
- **None**: No BLAST hits at all.

**Confidence score**: Proportion of hits (within the tier used) that agree with the majority-vote subtype. A confidence of 1.0 means all hits agree; 0.5 means half disagree.

**Third-pass results**:
- 57,914 total reads
- 57,912 with BLAST hits (99.997%)
- 54,253 high-confidence tier (>=98% identity, non-terminal)
- 3,659 low-confidence tier (<98% identity only)
- 2 reads with no BLAST hits

**Outputs**:
- `results/blast_validated.tsv` — **primary output**, one row per read
- `results/blast_all_hits.tsv` — all individual BLAST hits with parsed segment/subtype

## Third-Pass Results Summary

### Subtype Distribution

| Subtype | Reads | Notes |
|---------|-------|-------|
| H3N2 | ~47,000 | Dominant — consistent with seasonal circulation |
| H1N1 | ~6,800 | Second most common |
| H3N8 | ~1,500 | Equine/avian — detected across multiple sites |
| H5N1 | ~800 | Avian influenza — important for surveillance |
| H13N6 | ~400 | Gull-associated |
| H16N3 | ~300 | Shorebird-associated |
| H1N2 | ~200 | Swine-origin reassortant |
| H7N9 | ~50 | High-consequence avian subtype |
| H6N1 | ~30 | Avian |
| H10N5 | ~20 | Avian |
| Others | ~100+ | H5, H9N2, H4N6, H11N9, etc. |

### Segment Distribution

All 8 IAV segments represented. HA and NA segments are most informative for subtype assignment; internal segments (PB2, PB1, PA, NP, M, NS) reassort independently.

## Composition Visualization

After running the pipeline, generate subtype composition plots:

```bash
pip install matplotlib   # if not already installed
python plot_composition.py
```

Generates three figures in `results/`:

1. **`subtype_composition_by_site.png`** — Faceted stacked bars (3-column layout, one panel per site). X-axis = delivery dates, Y-axis = subtype proportion with total read count annotated above each bar.

2. **`subtype_composition_aggregate.png`** — Two-panel figure. Top: stacked subtype proportion bars across all sites combined. Bottom: total read count per delivery. Legend outside plot area to avoid overlap.

3. **`subtype_heatmap.png`** (300 DPI) — Sites x deliveries grid. Cell color = dominant subtype. White text showing `n=total` plus top 3 subtypes by percentage for each cell. Empty cells marked with `-`.

Major subtypes shown individually: H3N2, H1N1, H3N8, H5N1, H13N6, H16N3, H1N2. All others grouped as "Other".

## Multi-Pass BAM Support

The pipeline supports running on different alignment passes by changing the `bam_suffix` and `results_dir` config settings. Two config files are provided:

| Config | BAM suffix | Results directory | Description |
|--------|-----------|-------------------|-------------|
| `config.yaml` | `.third.filt.sorted.bam` | `results/` | Third-pass (tertiary) alignment — default |
| `config_first_pass.yaml` | `.first.filt.sorted.bam` | `results_first_pass/` | First-pass (initial) alignment |

```bash
# Third-pass run (default)
python run_pipeline.py

# First-pass run
python run_pipeline.py --config config_first_pass.yaml
```

The BLAST database (`data/blastdb/`) is shared between runs since it is independent of alignment pass.

**Note**: As of March 2025, first-pass BAM files (`.first.filt.sorted.bam`) are not uploaded to `gs://wwscanseq/viralmetagenomics/reports/`. Only third-pass BAMs are available in the GCS reports bucket. Running with `config_first_pass.yaml` will report "BAM not found" for all samples until first-pass BAMs are made available.

## Output Schema

### `blast_validated.tsv` (one row per read)

| Column | Description |
|--------|-------------|
| `read_id` | Original sequencer read name |
| `sample_ID` | Sample identifier (e.g., `ww035881`) |
| `City` | Sampling site city |
| `Date` | Sample collection date |
| `delivery_date` | Sequencing delivery batch date |
| `subspecies` | EsViritu classification (e.g., `t__Influenza A H3N2` or `t__H3N2`) |
| `accession` | EsViritu reference accession |
| `assigned_segment` | BLAST-assigned segment: HA, NA, PB2, PB1, PA, NP, M, or NS |
| `segment_confidence` | Proportion of hits agreeing on segment (0-1) |
| `assigned_subtype` | BLAST-assigned subtype (e.g., H3N2, H1N1) |
| `subtype_confidence` | Proportion of hits agreeing on subtype (0-1) |
| `identity_tier` | `high` (>=98% preferred hits), `low` (sub-threshold only), or `none` |
| `n_blast_hits` | Total BLAST hits for this read |
| `n_preferred_hits` | Hits passing both identity and terminal filters |
| `top_blast_hit` | Subject title of best BLAST hit (by bitscore) |
| `top_pident` | Percent identity of best hit |
| `top_evalue` | E-value of best hit |
| `top_bitscore` | Bitscore of best hit |
| `query_length` | Read length (bp) |
| `read_sequence` | Full nucleotide sequence |

### `blast_all_hits.tsv` (many rows per read)

Contains all columns from `blast_validated.tsv` metadata plus per-hit BLAST fields (`sseqid`, `stitle`, `pident`, `length`, `sstart`, `send`, `slen`, `evalue`, `bitscore`) and parsed `blast_segment`, `blast_subtype`, and flag columns (`above_identity_threshold`, `in_terminal_region`, `preferred_hit`).

## Configuration

All parameters are in `config.yaml`:

```yaml
# Paths
precomputed_virus_data: "../shiny_dashboard/dashboard_data/precomputed_virus.tsv.gz"
gcs_service_account: "../shiny_dashboard/gcs_service_account.json"
gcs_bucket: "wwscanseq"
gcs_prefix: "viralmetagenomics/reports"

# BAM file settings
bam_suffix: ".third.filt.sorted.bam"   # ".first.filt.sorted.bam" for first-pass
results_dir: "results"                   # override per run (e.g., "results_first_pass")

# Influenza A matching criteria
species_filter: "s__Alphainfluenzavirus influenzae"   # primary: filter by species column
subspecies_patterns:                                    # fallback: regex on subspecies column
  - "^t__Influenza A"
  - "^t__H\\d+N\\d+"

# BLAST database
blast_db_name: "influenza_a_complete"
entrez_query: '"Influenza A virus"[Organism] AND "complete cds"[Title]'
entrez_email: "azulli@stanford.edu"
entrez_batch_size: 10000

# BLAST parameters
max_target_seqs: 100
evalue: 1e-5
num_threads: 4

# Post-BLAST filtering (flags, not hard filters)
min_percent_identity: 98.0
terminal_filter_bp: 150
```

### Key Configuration Notes

- **`species_filter`** (added after initial run): Uses the species column as primary hit filter. This is more reliable than subspecies regex because EsViritu uses two naming conventions (`t__Influenza A H3N2` vs bare `t__H3N2`) and bare `t__H1`/`t__H2`/`t__H3` (without N) are actually Rotavirus subspecies. The species column (`s__Alphainfluenzavirus influenzae`) is consistent across all IAV entries.

- **`bam_suffix`**: Controls which BAM pass to extract reads from. Override in alternate config files for multi-pass analysis.

- **`results_dir`**: Controls where pipeline outputs are written. Use different directories per run to keep results separated.

- **`min_percent_identity`**: The 98% threshold flags hits but does not remove them. All BLAST hits are retained regardless of identity percentage.

## Requirements

- Python 3.9+
- BLAST+ (`blastn`, `makeblastdb`) on PATH
- GCS service account key at `../shiny_dashboard/gcs_service_account.json`
- Precomputed virus data at `../shiny_dashboard/dashboard_data/precomputed_virus.tsv.gz`

Python packages: `pandas`, `pysam`, `google-cloud-storage`, `biopython`, `pyyaml`

For composition plots: `matplotlib`, `numpy`

## Segment Reference

IAV has 8 genome segments. Subtype is defined by segments 4 (HA) and 6 (NA):

| Segment | Gene | Length (bp) | Subtype-defining? |
|---------|------|-------------|-------------------|
| 1 | PB2 | ~2,341 | No |
| 2 | PB1 | ~2,341 | No |
| 3 | PA | ~2,233 | No |
| 4 | **HA** | ~1,778 | **Yes (Hx)** |
| 5 | NP | ~1,565 | No |
| 6 | **NA** | ~1,413 | **Yes (Nx)** |
| 7 | M | ~1,027 | No |
| 8 | NS | ~890 | No |

Internal segments (PB2, PB1, PA, NP, M, NS) reassort independently of HA/NA. Subtype assignments from reads hitting internal segments reflect the closest BLAST match but are inherently less informative for subtyping than HA/NA reads.

## Directory Structure

```
influenza-blast-validation/
├── config.yaml                  # Pipeline parameters (third-pass default)
├── config_first_pass.yaml       # Alternate config for first-pass BAMs
├── requirements.txt             # Python dependencies
├── run_pipeline.py              # Orchestrator (--step and --config flags)
├── plot_composition.py          # Subtype composition visualization
├── scripts/
│   ├── __init__.py
│   ├── identify_hits.py         # Step 1: filter precomputed data
│   ├── extract_reads.py         # Step 2: BAM extraction from GCS
│   ├── build_blastdb.py         # Step 3: NCBI download + makeblastdb
│   ├── blast_validate.py        # Step 4: blastn + flag hits
│   └── assign_subtypes.py       # Step 5: segment/subtype assignment
├── data/
│   └── blastdb/                 # BLAST database files (~4 GB, gitignored)
├── results/                     # Third-pass pipeline outputs
│   ├── influenza_hits.tsv               # Step 1 output
│   ├── influenza_reads.fasta            # Step 2 output (gitignored)
│   ├── read_metadata.tsv                # Step 2 output (gitignored)
│   ├── blast_results_raw.tsv            # Step 4 raw output (gitignored)
│   ├── blast_results_flagged.tsv        # Step 4 flagged output (gitignored)
│   ├── blast_merged.tsv                 # Step 4 merged output (gitignored)
│   ├── blast_all_hits.tsv               # Step 5 all hits (gitignored)
│   ├── blast_validated.tsv              # Step 5 final per-read output
│   ├── subtype_composition_by_site.png  # Composition plots
│   ├── subtype_composition_aggregate.png
│   └── subtype_heatmap.png
└── results_first_pass/          # First-pass pipeline outputs (same structure)
```
