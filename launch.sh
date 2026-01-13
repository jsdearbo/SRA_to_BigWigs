#!/usr/bin/env bash
# Entry point orchestrating the RNA-seq pipeline.
# public data -> fastq -> trim -> align -> bigwig 

# Fail on first error, undefined variable, or pipeline failure for safer automation
set -euo pipefail

# Resolve the directory containing this script so relative paths stay stable
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
	cat <<'EOF'
Usage: launch.sh [options] [-- PROCESS_ARGS...]

Options:
  -c, --config PATH        Path to pipeline.env configuration file
  -m, --metadata PATH      Metadata CSV (overrides config for both stages)
	  --skip-download      Skip the download stage
	  --skip-process       Skip the processing stage
  -h, --help               Display this help message

Arguments after "--" are forwarded to processing_pipeline.sh.
EOF
}

# Track optional overrides and flags passed on the command line
CONFIG_OVERRIDE=""
METADATA_OVERRIDE=""
SKIP_DOWNLOAD=0
SKIP_PROCESS=0
declare -a PROCESS_ARGS=()

# Parse CLI flags for config/metadata overrides plus stage skip options
while (($#)); do
	case "${1}" in
		-c|--config)
			[[ $# -ge 2 ]] || { usage; exit 1; }
			CONFIG_OVERRIDE=${2}
			shift 2
			;;
		-m|--metadata)
			[[ $# -ge 2 ]] || { usage; exit 1; }
			METADATA_OVERRIDE=${2}
			shift 2
			;;
		--skip-download)
			SKIP_DOWNLOAD=1
			shift
			;;
		--skip-process)
			SKIP_PROCESS=1
			shift
			;;
		-h|--help)
			usage
			exit 0
			;;
		--)
			shift
			PROCESS_ARGS=("${@}")
			break
			;;
		*)
			usage
			exit 1
			;;
	esac
done

# Prefer an explicit --config path; otherwise fall back to repo-local config if present
CONFIG_ARG=()
if [[ -n ${CONFIG_OVERRIDE} ]]; then
	CONFIG_ARG=(--config "${CONFIG_OVERRIDE}")
else
	DEFAULT_CONFIG="${SCRIPT_DIR}/config/pipeline.env"
	if [[ -f ${DEFAULT_CONFIG} ]]; then
		CONFIG_ARG=(--config "${DEFAULT_CONFIG}")
	fi
fi

# Forward an optional metadata override to both downstream scripts
METADATA_ARG=()
if [[ -n ${METADATA_OVERRIDE} ]]; then
	METADATA_ARG=(--metadata "${METADATA_OVERRIDE}")
fi

# Trigger the download stage unless explicitly skipped
if (( ! SKIP_DOWNLOAD )); then
	echo "Step 1: Downloading missing SRA files..."
	"${SCRIPT_DIR}/download_data.sh" "${CONFIG_ARG[@]}" "${METADATA_ARG[@]}"
fi

# Launch processing with any extra arguments passed after "--" unless skipped
if (( ! SKIP_PROCESS )); then
	echo "Step 2: Processing datasets..."
	"${SCRIPT_DIR}/processing_pipeline.sh" "${CONFIG_ARG[@]}" "${METADATA_ARG[@]}" "${PROCESS_ARGS[@]}"
fi

echo "Pipeline finished successfully!"