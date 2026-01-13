#!/usr/bin/env bash
set -euo pipefail

# Determine script location and load shared helpers for logging/config.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# Render CLI help text for end users.
usage() {
    cat <<'EOF'
Usage: download_data.sh [options]

Options:
  -c, --config PATH     Path to pipeline.env configuration file
  -m, --metadata PATH   Metadata CSV describing runs to fetch
  -o, --output  PATH    Directory where SRA files will be stored
  -n, --dry-run         Show planned downloads without fetching
  -h, --help            Display this help message

Values fall back to config defaults, then repository defaults.
EOF
}

# Capture command-line overrides before loading configs.
CONFIG_FILE=""
METADATA_OVERRIDE=""
OUTPUT_OVERRIDE=""
DRY_RUN=0

# Parse flags for config, metadata, output, and dry-run mode.
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
        -o|--output)
            [[ $# -ge 2 ]] || die "Missing value for ${1}"
            OUTPUT_OVERRIDE=${2}
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=1
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

# Load configuration values either from CLI or the default location.
if [[ -n ${CONFIG_FILE} ]]; then
    load_config "${CONFIG_FILE}"
else
    load_config "$CONFIG_FILE_DEFAULT"
fi

# Resolve defaults and overrides for metadata and output directories.
DEFAULT_METADATA="${PIPELINE_ROOT}/metadata.csv"
DEFAULT_OUTPUT="${PIPELINE_ROOT}/sra_files"

METADATA_INPUT="${METADATA_OVERRIDE:-${METADATA_FILE:-$DEFAULT_METADATA}}"
OUTPUT_INPUT="${OUTPUT_OVERRIDE:-${SRA_OUTPUT_DIR:-$DEFAULT_OUTPUT}}"

METADATA_FILE=$(resolve_path "${METADATA_INPUT}")
OUTPUT_DIR=$(resolve_path "${OUTPUT_INPUT}")

# Ensure prerequisites and filesystem layout are ready.
require_tools prefetch python3

[[ -f ${METADATA_FILE} ]] || die "Metadata file not found: ${METADATA_FILE}"
ensure_directory "${OUTPUT_DIR}"

info "Using metadata ${METADATA_FILE}"
info "Writing SRA files to ${OUTPUT_DIR}"

# Extract required metadata columns via Python for robust CSV parsing.
read_metadata() {
    python3 - "$METADATA_FILE" <<'PY'
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
        cell = (row.get("CellType") or "").strip()
        run = (row.get("RunID") or "").strip()
        desc = (row.get("Description") or "").replace("\n", " ").strip()
        if not run:
            continue
        print("\t".join([cell, run, desc]))
PY
}

# Iterate over metadata rows, downloading only missing SRA files.
while IFS=$'\t' read -r cell_type run_id description; do
    [[ -n ${run_id} ]] || continue

    info "Target: ${cell_type:-unknown} | ID: ${run_id}"
    [[ -z ${description} ]] || info "Desc: ${description}"

    target_dir="${OUTPUT_DIR}/${run_id}"
    target_file="${target_dir}/${run_id}.sra}"
    alt_target="${OUTPUT_DIR}/${run_id}.sra"

    if [[ -f ${target_file} || -f ${alt_target} ]]; then
        info "[SKIP] ${run_id} already present"
        continue
    fi

    if ((DRY_RUN)); then
        info "[DRY-RUN] Would fetch ${run_id}"
        continue
    fi

    info "[DOWNLOADING] Fetching ${run_id}"
    prefetch --max-size 50G -O "${OUTPUT_DIR}" "${run_id}"
done < <(read_metadata)

# Summarize the stage once all entries have been processed.
info "All downloads processed based on ${METADATA_FILE}"