# Influenza A BLAST Validation Pipeline

Extracts Influenza A virus-classified reads from EsViritu wastewater metagenomics BAM files, validates them via BLAST against all complete IAV coding sequences in NCBI (~1.33M sequences), assigns segment identity and subtype by majority vote across top BLAST hits, and reports per-read confidence scores.

## Quick Start

```bash
cd influenza-blast-validation
pip install -r requirements.txt    # pandas, pysam, google-cloud-storage, biopython, pyyaml
python run_pipeline.py             # runs all 5 steps
```

Run individual steps:

```bash
python run_pipeline.py --step identify   # step 1 only
python run_pipeline.py --step blastdb    # step 3 only (can run in parallel with step 2)
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

Reads `precomputed_virus.tsv.gz` from the metagenomics dashboard and filters for EsViritu subspecies matching `^t__Influenza A` (regex). Captures all IAV subtypes: H3N2, H1N1, H3N8, H5N1, H13N2, etc.

**Output**: `results/influenza_hits.tsv`

### Step 2: Extract Reads

For each hit, downloads the corresponding BAM file from GCS (`gs://wwscanseq/viralmetagenomics/reports/`), extracts aligned reads by accession number, and writes sequences + metadata.

**Outputs**:
- `results/influenza_reads.fasta` — query sequences for BLAST
- `results/read_metadata.tsv` — per-read metadata (sample, site, date, EsViritu classification)

### Step 3: Build BLAST Database

Downloads all Influenza A virus complete CDS sequences from NCBI Entrez in batches of 10,000 using the search history API. Builds a local BLAST nucleotide database.

- **Query**: `"Influenza A virus"[Organism] AND "complete cds"[Title]`
- **Database size**: ~1.33M sequences, ~2.4 GB FASTA, ~1.5–2 GB BLAST DB
- Skips download if FASTA already exists (>1 GB); skips DB build if DB files already exist
- Delete `data/blastdb/genomes.fasta` or DB files to force rebuild

**Output**: `data/blastdb/influenza_a_complete.*`

### Step 4: BLAST Validation

Runs `blastn` with up to 100 target sequences per read. All hits are retained and flagged (not filtered):

| Flag | Column | Meaning |
|------|--------|---------|
| Identity threshold | `above_identity_threshold` | True if percent identity >= 98% |
| Terminal region | `in_terminal_region` | True if alignment falls entirely within the first or last 150 bp of the subject sequence |
| Preferred hit | `preferred_hit` | True if above threshold AND not in terminal region |

The 150 bp terminal filter addresses conserved UTR/signal peptide/transmembrane regions at segment ends that can produce cross-subtype BLAST matches (IAV UTRs are 20–46 bp, with an additional ~100 bp of moderately conserved coding sequence beyond them).

Skips the BLAST search if `blast_results_raw.tsv` already exists. Delete it to force re-run.

**Outputs**:
- `results/blast_results_raw.tsv` — raw BLAST output (tab-separated, no header)
- `results/blast_results_flagged.tsv` — all hits with flag columns
- `results/blast_merged.tsv` — hits left-joined to read metadata (many rows per read)

### Step 5: Assign Subtypes

Parses segment and subtype from each BLAST hit's subject title, then assigns per-read values by majority vote.

**Segment identification** — parsed from BLAST subject titles using (in priority order):
1. Segment number: `"segment 4"` → HA
2. Gene abbreviation: `"(HA)"`, `"(NA)"`, `"(PB2)"`
3. Full gene name: `"hemagglutinin"`, `"neuraminidase"`

**Subtype identification** — parsed from strain name parenthetical:
- `A/swine/Minnesota/.../2026(H1N2)` → H1N2
- `A/duck/.../2024(H5)` → H5 (NA not annotated)
- Mixed subtypes are excluded

**Tiered assignment**:
- **High tier**: read has preferred hits (≥98% identity, non-terminal). Majority vote uses only preferred hits.
- **Low tier**: read has BLAST hits but none are preferred. Majority vote uses all hits. These reads are still assigned a subtype but should be interpreted with caution.
- **None**: no BLAST hits at all.

**Confidence score**: proportion of hits (within the tier used) that agree with the majority-vote subtype. A confidence of 1.0 means all hits agree; 0.5 means half disagree.

**Outputs**:
- `results/blast_validated.tsv` — **primary output**, one row per read
- `results/blast_all_hits.tsv` — all individual BLAST hits with parsed segment/subtype

## Output Schema

### `blast_validated.tsv` (one row per read)

| Column | Description |
|--------|-------------|
| `read_id` | Original sequencer read name |
| `sample_ID` | Sample identifier (e.g., `ww035881`) |
| `City` | Sampling site city |
| `Date` | Sample collection date |
| `delivery_date` | Sequencing delivery batch date |
| `subspecies` | EsViritu classification (e.g., `t__Influenza A H3N2`) |
| `accession` | EsViritu reference accession |
| `assigned_segment` | BLAST-assigned segment: HA, NA, PB2, PB1, PA, NP, M, or NS |
| `segment_confidence` | Proportion of hits agreeing on segment (0–1) |
| `assigned_subtype` | BLAST-assigned subtype (e.g., H3N2, H1N1) |
| `subtype_confidence` | Proportion of hits agreeing on subtype (0–1) |
| `identity_tier` | `high` (≥98% preferred hits), `low` (sub-threshold only), or `none` |
| `n_blast_hits` | Total BLAST hits for this read |
| `n_preferred_hits` | Hits passing both identity and terminal filters |
| `top_blast_hit` | Subject title of best BLAST hit (by bitscore) |
| `top_pident` | Percent identity of best hit |
| `top_evalue` | E-value of best hit |
| `top_bitscore` | Bitscore of best hit |
| `query_length` | Read length (bp) |
| `read_sequence` | Full nucleotide sequence |

### `blast_all_hits.tsv` (many rows per read)

Contains all columns from `blast_validated.tsv` metadata plus per-hit BLAST fields (`sseqid`, `stitle`, `pident`, `length`, `sstart`, `send`, `slen`, `evalue`, `bitscore`) and parsed `blast_segment`, `blast_subtype`, and flag columns.

## Configuration

All parameters are in `config.yaml`:

```yaml
# Subspecies patterns (regex)
subspecies_patterns:
  - "^t__Influenza A"

# BLAST database
entrez_query: '"Influenza A virus"[Organism] AND "complete cds"[Title]'
entrez_batch_size: 10000

# BLAST parameters
max_target_seqs: 100
evalue: 1e-5
num_threads: 4

# Post-BLAST filtering (flags, not hard filters)
min_percent_identity: 98.0
terminal_filter_bp: 150
```

## Requirements

- Python 3.9+
- BLAST+ (`blastn`, `makeblastdb`) on PATH
- GCS service account key at `../shiny_dashboard/gcs_service_account.json`
- Precomputed virus data at `../shiny_dashboard/dashboard_data/precomputed_virus.tsv.gz`

Python packages: `pandas`, `pysam`, `google-cloud-storage`, `biopython`, `pyyaml`

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
├── config.yaml              # All pipeline parameters
├── requirements.txt         # Python dependencies
├── run_pipeline.py          # Orchestrator (--step flag for individual steps)
├── scripts/
│   ├── identify_hits.py     # Step 1: filter precomputed data
│   ├── extract_reads.py     # Step 2: BAM extraction from GCS
│   ├── build_blastdb.py     # Step 3: NCBI download + makeblastdb
│   ├── blast_validate.py    # Step 4: blastn + flag hits
│   └── assign_subtypes.py   # Step 5: segment/subtype assignment
├── data/
│   └── blastdb/             # BLAST database files (~4 GB)
└── results/
    ├── influenza_hits.tsv           # Step 1 output
    ├── influenza_reads.fasta        # Step 2 output
    ├── read_metadata.tsv            # Step 2 output
    ├── blast_results_raw.tsv        # Step 4 raw output
    ├── blast_results_flagged.tsv    # Step 4 flagged output
    ├── blast_merged.tsv             # Step 4 merged output
    ├── blast_all_hits.tsv           # Step 5 all hits with parsed fields
    └── blast_validated.tsv          # Step 5 final per-read output
```
