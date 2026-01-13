# RNA-Seq Processing Pipeline

Pipeline for downloading public SRA accessions, generating cleaned FASTQ files, aligning with STAR, and producing normalized BigWig coverage tracks. The workflow is driven by a metadata CSV and configurable through a simple environment file so it can be shared with lab members or versioned on GitHub.

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

Install these tools in a shared conda environment (recommended) or expose their executables on `PATH`. Absolute paths can be specified in the config file when needed.

## Repository Layout
- `download_data.sh` – fetches `.sra` files listed in the metadata
- `processing_pipeline.sh` – converts SRA to FASTQ, trims, aligns, and emits BAM/BigWig
- `launch.sh` – convenience wrapper that runs download + processing back-to-back
- `metadata.csv` – example metadata (CellType, RunID, StudyID, Description, SeqType, Year)
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

## Version Control Tips
- Commit `config/pipeline.env.example`, but keep your personalized `config/pipeline.env` out of version control (add it to `.gitignore`).
- Document environment creation (e.g., `conda env export --from-history > environment.yml`) so lab members can reproduce dependencies.
- Use branches or pull requests to propose changes to parameters or tool versions and capture the history of updates.
