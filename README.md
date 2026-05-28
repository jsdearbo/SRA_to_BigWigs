# RNA-Seq Processing Pipeline

Pipeline for downloading public SRA accessions, generating cleaned FASTQ files, aligning with STAR, and producing normalized BigWig coverage tracks. The workflow is driven by a metadata CSV and configurable through a simple environment file so it can be shared with lab members or versioned on GitHub.

## Overview

```text
    launch.sh
    │
    ├─ Step 1  download_data.sh
    │            reads metadata.csv
    │            └─ prefetch ──► sra_files/{RunID}/{RunID}.sra
    │
    └─ Step 2  processing_pipeline.sh  (one pass per RunID in metadata.csv)
                 │
                 ├─ fasterq-dump ──► fastq/{srr}_1.fastq
                 │                   fastq/{srr}_2.fastq
                 │
                 ├─ Trimmomatic  ──► trimmed/{srr}_1_trimmed.fastq
                 │                   trimmed/{srr}_2_trimmed.fastq
                 │
                 ├─ STAR         ──► bam/{srr}.bam  bam/{srr}.bam.bai
                 │                   logs/
                 │
                 └─ bamCoverage  ──► bw/{srr}.bw
```

## Features
- Single entrypoint with optional stage skipping for incremental runs
- Configurable paths, tool locations, and resource usage via `config/pipeline.env`
- Structured outputs grouped by cell type with logs, intermediates, and tracks
- Safe defaults with dry-run support for SRA downloads and per-sample filtering

## Requirements
- Bash 4+
- SRA Toolkit (`prefetch`, `fasterq-dump`)
- STAR aligner
- Trimmomatic
- samtools
- deepTools (`bamCoverage`)
- Python 3.7+

Install via the provided `environment.yml`:
```bash
conda env create -f environment.yml
conda activate rnaseq_tools
```
Or expose the executables on `PATH` directly. Absolute paths can be specified in the config file when needed.

## Repository Layout
- `download_data.sh` – fetches `.sra` files listed in the metadata
- `processing_pipeline.sh` – converts SRA to FASTQ, trims, aligns, and emits BAM/BigWig
- `launch.sh` – convenience wrapper that runs download + processing back-to-back
- `metadata.csv` – example metadata (CellType, RunID, StudyID, Description, SeqType, Year)
- `environment.yml` – conda environment definition for all required tools
- `config/pipeline.env.example` – sample configuration to copy and customize
- `lib/common.sh` – shared helper functions for logging and path handling

## Setup
1. Clone or download this directory into your workspace.
2. Copy the sample config and edit it for your environment:
	```bash
	cp config/pipeline.env.example config/pipeline.env
	$EDITOR config/pipeline.env
	```
	Update tool paths (STAR, Trimmomatic, etc.), adjust thread counts, and set output locations relative to the repo or as absolute paths.
3. Ensure your metadata CSV follows the schema shown in `metadata.csv`.

## Usage
Run the full pipeline:
```bash
./launch.sh
```
Key options:
- `./launch.sh --config config/pipeline.env` – use a specific configuration file
- `./launch.sh --metadata metadata.csv` – override the metadata file for both stages
- `./launch.sh --skip-download` – reuse existing `.sra` files
- `./launch.sh --skip-process` – only refresh downloads
- `./launch.sh -- --sample SRR123456` – forward additional arguments to `processing_pipeline.sh`

Stage entrypoints expose more granular controls:
- Dry-run SRA downloads: `./download_data.sh --dry-run`
- Limit processing to selected runs: `./processing_pipeline.sh --sample SRR17143399`
- Skip individual steps (e.g., `--skip-trim`) when intermediates already exist
- Process single-end libraries: `./processing_pipeline.sh --single-end`

## Outputs
For each `CellType` row in the metadata, the pipeline creates a directory under `BASE_DIR` containing:
- `fastq/` – raw FASTQ files
- `trimmed/` – adapter-trimmed FASTQs
- `bam/` – coordinate-sorted BAMs plus indexes
- `bw/` – CPM-normalized BigWig tracks
- `logs/` – STAR logs and splice junction summaries
- `tmp/` – scratch space (removed automatically unless `KEEP_TMP=1`)

## Metadata Format
`metadata.csv` must include at least the following columns:

| CellType | RunID | Description | ... |
|----------|-------|-------------|-----|

Additional columns are ignored by the scripts, so you can extend the file with project-specific annotations.

## Troubleshooting

**Resuming a partial run:** The pipeline is idempotent — all steps check for existing outputs before running. Simply re-run the same command to pick up where it left off.

**Isolating a single sample:** Use `--sample SRR######` with `processing_pipeline.sh` to process one run at a time. This is useful for debugging or re-running a sample that failed.

**Inspecting failures:** Per-sample logs are written to `results/<CellType>/logs/`. STAR alignment logs (`Log.final.out`) and Trimmomatic output are the first places to look.

**Missing or stale FASTQs:** Pass `--skip-fastq` or `--skip-trim` when intermediates already exist but you want to re-run downstream steps only.

## Testing

A smoke test suite is included that validates config loading, metadata parsing, and CLI flags without requiring any bioinformatics tools:

```bash
bash tests/test_pipeline.sh
```
