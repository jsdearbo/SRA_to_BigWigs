#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# Print CLI usage information for end users
usage() {
    cat <<'EOF'
Usage: processing_pipeline.sh [options]

Options:
  -c, --config PATH           Path to pipeline.env configuration file
  -m, --metadata PATH         Metadata CSV describing samples
  -b, --base-dir PATH         Output directory root (per cell type subfolders)
      --sra-dir PATH          Directory containing SRA files
      --star-index PATH       STAR genome index directory
      --star-exec PATH        STAR executable path or name
      --trimmomatic PATH      Trimmomatic executable path or name
      --adapters PATH         Adapter FASTA for Trimmomatic
  -t, --threads INT           Alignment threads (default 6)
  -i, --io-threads INT        I/O threads for fasterq-dump and bamCoverage (default 4)
      --sample RUN_ID         Restrict processing to specific run ID (repeatable)
      --skip-fastq            Skip FASTQ generation step (requires existing outputs)
      --skip-trim             Skip trimming step
      --skip-align            Skip STAR alignment
      --skip-bw               Skip BigWig generation
      --keep-tmp              Retain per-sample temporary directories
  -h, --help                  Display this help message
EOF
}

# Collect CLI overrides and bookkeeping flags
CONFIG_FILE=""
METADATA_OVERRIDE=""
BASE_OVERRIDE=""
SRA_OVERRIDE=""
STAR_INDEX_OVERRIDE=""
STAR_EXEC_OVERRIDE=""
TRIMMOMATIC_OVERRIDE=""
ADAPTERS_OVERRIDE=""
THREADS_OVERRIDE=""
IO_THREADS_OVERRIDE=""
KEEP_TMP_OVERRIDE=""
declare -a SAMPLE_FILTER=()

SKIP_FASTQ=0
SKIP_TRIM=0
SKIP_ALIGN=0
SKIP_BW=0

# Parse command-line arguments to populate overrides and skip flags
while (($#)); do
    case "${1}" in
        -c|--config)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            CONFIG_FILE=$(resolve_path "${2}")
            shift 2
            ;;
        -m|--metadata)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            METADATA_OVERRIDE=${2}
            shift 2
            ;;
        -b|--base-dir)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            BASE_OVERRIDE=${2}
            shift 2
            ;;
        --sra-dir)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            SRA_OVERRIDE=${2}
            shift 2
            ;;
        --star-index)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            STAR_INDEX_OVERRIDE=${2}
            shift 2
            ;;
        --star-exec)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            STAR_EXEC_OVERRIDE=${2}
            shift 2
            ;;
        --trimmomatic)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            TRIMMOMATIC_OVERRIDE=${2}
            shift 2
            ;;
        --adapters)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            ADAPTERS_OVERRIDE=${2}
            shift 2
            ;;
        -t|--threads)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            THREADS_OVERRIDE=${2}
            shift 2
            ;;
        -i|--io-threads)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            IO_THREADS_OVERRIDE=${2}
            shift 2
            ;;
        --sample)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            SAMPLE_FILTER+=("${2}")
            shift 2
            ;;
        --skip-fastq)
            SKIP_FASTQ=1
            shift
            ;;
        --skip-trim)
            SKIP_TRIM=1
            shift
            ;;
        --skip-align)
            SKIP_ALIGN=1
            shift
            ;;
        --skip-bw)
            SKIP_BW=1
            shift
            ;;
        --keep-tmp)
            KEEP_TMP_OVERRIDE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage
            die "Unknown argument: ${1}"
            ;;
    esac
done

# Load configuration from an explicit file or fall back to the default
if [[ -n ${CONFIG_FILE} ]]; then
    load_config "${CONFIG_FILE}"
else
    load_config "$CONFIG_FILE_DEFAULT"
fi

# Establish default values used when config or CLI omit a setting
DEFAULT_METADATA="${PIPELINE_ROOT}/metadata.csv"
DEFAULT_BASE_DIR="${PIPELINE_ROOT}"
DEFAULT_SRA_DIR="${PIPELINE_ROOT}/sra_files"
DEFAULT_STAR_EXEC="/mnt/netfiles/conda_envs/rna-seq/bin/STAR"
DEFAULT_STAR_INDEX="/mnt/netfiles/genomes/STAR/mm10/star_index"
DEFAULT_TRIMMOMATIC_EXEC="/mnt/netfiles/conda_envs/rna-seq/bin/trimmomatic"
DEFAULT_ADAPTERS="/mnt/netfiles/conda_envs/rna-seq/share/trimmomatic/adapters/TruSeq3-PE.fa"
DEFAULT_ALIGN_THREADS=6
DEFAULT_IO_THREADS=4
DEFAULT_STAR_OVERHANG=100
DEFAULT_EFFECTIVE_GENOME_SIZE=2652783500
DEFAULT_BW_BIN_SIZE=10
DEFAULT_BW_NORMALIZATION=CPM
DEFAULT_JAVA_OPTIONS="-Xmx8G"

# Resolve inputs by applying CLI overrides, config values, or defaults
METADATA_INPUT="${METADATA_OVERRIDE:-${METADATA_FILE:-$DEFAULT_METADATA}}"
BASE_INPUT="${BASE_OVERRIDE:-${BASE_DIR:-$DEFAULT_BASE_DIR}}"
SRA_INPUT="${SRA_OVERRIDE:-${SRA_SOURCE_DIR:-$DEFAULT_SRA_DIR}}"
STAR_INDEX_INPUT="${STAR_INDEX_OVERRIDE:-${STAR_INDEX:-$DEFAULT_STAR_INDEX}}"
STAR_EXEC_INPUT="${STAR_EXEC_OVERRIDE:-${STAR_EXEC:-$DEFAULT_STAR_EXEC}}"
TRIMMOMATIC_INPUT="${TRIMMOMATIC_OVERRIDE:-${TRIMMOMATIC_EXEC:-$DEFAULT_TRIMMOMATIC_EXEC}}"
ADAPTERS_INPUT="${ADAPTERS_OVERRIDE:-${ADAPTERS:-$DEFAULT_ADAPTERS}}"

ALIGN_THREADS_VAL="${THREADS_OVERRIDE:-${ALIGN_THREADS:-$DEFAULT_ALIGN_THREADS}}"
IO_THREADS_VAL="${IO_THREADS_OVERRIDE:-${IO_THREADS:-$DEFAULT_IO_THREADS}}"
STAR_OVERHANG_VAL="${STAR_OVERHANG:-$DEFAULT_STAR_OVERHANG}"
EFFECTIVE_GENOME_SIZE_VAL="${EFFECTIVE_GENOME_SIZE:-$DEFAULT_EFFECTIVE_GENOME_SIZE}"
BW_BIN_SIZE_VAL="${BW_BIN_SIZE:-$DEFAULT_BW_BIN_SIZE}"
BW_NORMALIZATION_VAL="${BW_NORMALIZATION:-$DEFAULT_BW_NORMALIZATION}"
JAVA_OPTIONS_VAL="${_JAVA_OPTIONS:-$DEFAULT_JAVA_OPTIONS}"

KEEP_TMP_VAL=${KEEP_TMP_OVERRIDE:-${KEEP_TMP:-0}}

# Convert user inputs into absolute paths and resolvable executables
METADATA_FILE=$(resolve_path "${METADATA_INPUT}")
BASE_DIR=$(resolve_path "${BASE_INPUT}")
SRA_SOURCE_DIR=$(resolve_path "${SRA_INPUT}")
STAR_INDEX=$(resolve_path "${STAR_INDEX_INPUT}")
ADAPTERS=$(resolve_path "${ADAPTERS_INPUT}")

STAR_EXEC=$(resolve_executable "${STAR_EXEC_INPUT}")
TRIMMOMATIC_EXEC=$(resolve_executable "${TRIMMOMATIC_INPUT}")
SAMTOOLS_BIN=$(resolve_executable "${SAMTOOLS_BIN:-samtools}")
BAMCOVERAGE_BIN=$(resolve_executable "${BAMCOVERAGE_BIN:-bamCoverage}")
FASTERQ_DUMP_BIN=$(resolve_executable "${FASTERQ_DUMP_BIN:-fasterq-dump}")

ALIGN_THREADS=${ALIGN_THREADS_VAL}
IO_THREADS=${IO_THREADS_VAL}
STAR_OVERHANG=${STAR_OVERHANG_VAL}
EFFECTIVE_GENOME_SIZE=${EFFECTIVE_GENOME_SIZE_VAL}
BW_BIN_SIZE=${BW_BIN_SIZE_VAL}
BW_NORMALIZATION=${BW_NORMALIZATION_VAL}
export _JAVA_OPTIONS="${JAVA_OPTIONS_VAL}"

# Confirm required tooling and inputs are ready before proceeding
require_tools python3

[[ -f ${METADATA_FILE} ]] || die "Metadata file not found: ${METADATA_FILE}"
[[ -d ${SRA_SOURCE_DIR} ]] || die "SRA directory not found: ${SRA_SOURCE_DIR}"
[[ -d ${STAR_INDEX} ]] || die "STAR index directory not found: ${STAR_INDEX}"
[[ -f ${ADAPTERS} ]] || die "Adapter FASTA not found: ${ADAPTERS}"

# Expand optional extra argument strings into arrays for downstream commands
IFS=' ' read -r -a STAR_EXTRA_ARGS_ARR <<< "${STAR_EXTRA_ARGS:-}"
IFS=' ' read -r -a BAMCOVERAGE_EXTRA_ARGS_ARR <<< "${BAMCOVERAGE_EXTRA_ARGS:-}"
IFS=' ' read -r -a TRIMMOMATIC_EXTRA_ARGS_ARR <<< "${TRIMMOMATIC_EXTRA_ARGS:-}"

# Create the top-level output directory if it does not exist
ensure_directory "${BASE_DIR}"

# Log key paths and configuration for debugging and reproducibility
info "Pipeline root: ${PIPELINE_ROOT}"
info "Metadata: ${METADATA_FILE}"
info "Output base: ${BASE_DIR}"
info "SRA source: ${SRA_SOURCE_DIR}"
info "STAR index: ${STAR_INDEX}"

# Decide whether to handle a sample based on optional run filters
should_process_run() {
    local run_id=$1
    if ((${#SAMPLE_FILTER[@]} == 0)); then
        return 0
    fi
    local candidate
    for candidate in "${SAMPLE_FILTER[@]}"; do
        if [[ ${candidate} == "${run_id}" ]]; then
            return 0
        fi
    done
    return 1
}

# Stream metadata rows (cell type, run, description) via Python parsing
metadata_rows() {
    python3 - "${METADATA_FILE}" <<'PY'
import csv
import sys

metadata_path = sys.argv[1]

with open(metadata_path, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    required = {"CellType", "RunID", "Description"}
    if missing := (required - set(reader.fieldnames or [])):
        sys.stderr.write(f"Missing columns in metadata: {', '.join(sorted(missing))}\n")
        sys.exit(2)
    for row in reader:
        run = (row.get("RunID") or "").strip()
        if not run:
            continue
        cell = (row.get("CellType") or "").strip()
        desc = (row.get("Description") or "").replace("\n", " ").strip()
        print("\t".join([cell, run, desc]))
PY
}

# Perform all processing steps for a single SRR sample (Main pipeline function)
process_sample() {
    local cell_type=$1
    local srr=$2
    local description=$3

    # Respect optional --sample filters before doing any work
    should_process_run "${srr}" || {
        info "Skipping ${srr} (not in filter)"
        return 0
    }

    # Create a per-cell-type directory layout for all outputs
    local out_root="${BASE_DIR}/${cell_type:-unspecified}"
    ensure_directory "${out_root}"
    ensure_directory "${out_root}/fastq"
    ensure_directory "${out_root}/trimmed"
    ensure_directory "${out_root}/bam"
    ensure_directory "${out_root}/bw"
    ensure_directory "${out_root}/logs"
    ensure_directory "${out_root}/tmp"

    # Prepare a temporary workspace used by fasterq-dump and STAR
    local sample_tmp="${out_root}/tmp/${srr}"
    ensure_directory "${sample_tmp}"
    ensure_directory "${sample_tmp}/fasterq"
    ensure_directory "${sample_tmp}/star"
    # Remove any stale STAR tmp directory so the aligner can recreate it
    rm -rf "${sample_tmp}/star/tmp"

    # Conditionally delete temporary directories on exit unless requested otherwise
    if (( KEEP_TMP_VAL )); then
        info "Temporary files retained at ${sample_tmp}"
    else
        trap 'rm -rf "${sample_tmp}"; trap - RETURN' RETURN
    fi

    # Emit a formatted header so logs are easy to scan
    info "=================================================="
    info "Processing: ${srr} (${cell_type:-NA})"
    [[ -z ${description} ]] || info "Desc: ${description}"
    info "=================================================="

    # Generate paired FASTQ files from the SRA archive if needed
    local fq1="${out_root}/fastq/${srr}_1.fastq"
    local fq2="${out_root}/fastq/${srr}_2.fastq"

    if (( SKIP_FASTQ )); then
        info "[FASTQ] Skip requested"
        [[ -f ${fq1} && -f ${fq2} ]] || die "FASTQ files missing for ${srr} while --skip-fastq is set"
    elif [[ ! -f ${fq1} || ! -f ${fq2} ]]; then
        info "[FASTQ] Dumping ${srr}"
        local sra_path="${SRA_SOURCE_DIR}/${srr}/${srr}.sra"
        [[ -f ${sra_path} ]] || sra_path="${SRA_SOURCE_DIR}/${srr}.sra"
        [[ -f ${sra_path} ]] || die "SRA file not found for ${srr}. Run download_data.sh first."

        "${FASTERQ_DUMP_BIN}" --split-3 \
            --threads "${IO_THREADS}" \
            --outdir "${out_root}/fastq" \
            --temp "${sample_tmp}/fasterq" \
            "${sra_path}"
        sync
    else
        info "[FASTQ] Exists — skipping"
    fi

    # Run Trimmomatic to clean adapters/low-quality bases if outputs are absent
    local trim_fq1="${out_root}/trimmed/${srr}_1_trimmed.fastq"
    local trim_fq2="${out_root}/trimmed/${srr}_2_trimmed.fastq"

    if (( SKIP_TRIM )); then
        info "[TRIM] Skip requested"
        [[ -f ${trim_fq1} && -f ${trim_fq2} ]] || die "Trimmed FASTQ files missing for ${srr} while --skip-trim is set"
    elif [[ ! -f ${trim_fq1} || ! -f ${trim_fq2} ]]; then
        info "[TRIM] Trimming reads"
        "${TRIMMOMATIC_EXEC}" PE -threads "${ALIGN_THREADS}" \
            "${fq1}" "${fq2}" \
            "${trim_fq1}" "${sample_tmp}/1_unpaired.fq" \
            "${trim_fq2}" "${sample_tmp}/2_unpaired.fq" \
            ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            "${TRIMMOMATIC_EXTRA_ARGS_ARR[@]}"
        rm -f "${sample_tmp}/"*_unpaired.fq 2>/dev/null || true
        sync
    else
        info "[TRIM] Exists — skipping"
    fi

    # Align trimmed reads with STAR and produce a coordinate-sorted BAM
    local final_bam="${out_root}/bam/${srr}.bam"

    if (( SKIP_ALIGN )); then
        info "[STAR] Skip requested"
        [[ -f ${final_bam} ]] || die "BAM missing for ${srr} while --skip-align is set"
    elif [[ ! -f ${final_bam} ]]; then
        info "[STAR] Aligning ${srr}"
        "${STAR_EXEC}" \
            --runThreadN "${ALIGN_THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${trim_fq1}" "${trim_fq2}" \
            --outFileNamePrefix "${sample_tmp}/star/" \
            --outTmpDir "${sample_tmp}/star/tmp" \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbOverhang "${STAR_OVERHANG}" \
            --outBAMcompression 1 \
            --outSAMattributes All \
            "${STAR_EXTRA_ARGS_ARR[@]}"

        mv "${sample_tmp}/star/Aligned.sortedByCoord.out.bam" "${final_bam}"
        mv "${sample_tmp}/star/Log."* "${out_root}/logs/" 2>/dev/null || true
        mv "${sample_tmp}/star/SJ.out.tab" "${out_root}/logs/" 2>/dev/null || true

        "${SAMTOOLS_BIN}" index "${final_bam}"
        sync
    else
        info "[STAR] BAM exists — skipping"
    fi

    # Convert BAM coverage to a BigWig track when requested
    local bw_out="${out_root}/bw/${srr}.bw"

    if (( SKIP_BW )); then
        info "[BW] Skip requested"
        [[ -f ${bw_out} ]] || warn "BigWig missing for ${srr} while --skip-bw is set"
    elif [[ ! -f ${bw_out} ]]; then
        info "[BW] Generating BigWig"
        "${BAMCOVERAGE_BIN}" \
            -b "${final_bam}" \
            -o "${bw_out}" \
            --binSize "${BW_BIN_SIZE}" \
            --normalizeUsing "${BW_NORMALIZATION}" \
            --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" \
            --numberOfProcessors "${IO_THREADS}" \
            "${BAMCOVERAGE_EXTRA_ARGS_ARR[@]}"
        sync
    else
        info "[BW] Exists — skipping"
    fi

    # Clean up temporary directories unless the user asked to keep them
    if (( KEEP_TMP_VAL == 0 )); then
        rm -rf "${sample_tmp}"
    fi

    # Mark this run as successfully processed
    info "[DONE] ${srr}"
}

# Iterate over every metadata entry and process matching samples
metadata_rows | while IFS=$'\t' read -r cell_type run_id description; do
    process_sample "${cell_type}" "${run_id}" "${description}"
done

# Report successful completion of all requested tasks
info "All pipeline tasks complete."