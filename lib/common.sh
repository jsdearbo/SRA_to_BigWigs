#!/usr/bin/env bash
# Shared helpers for pipeline scripts.

# Prevent this helper from being sourced multiple times.
if [[ -n ${PIPELINE_COMMON_SH:-} ]]; then
    return 0
fi
export PIPELINE_COMMON_SH=1

# Determine the absolute directory that contains this file, following symlinks.
_current_script_dir() {
    local src
    src=${BASH_SOURCE[0]}
    while [[ -h "$src" ]]; do
        local dir
        dir=$(cd -P "$(dirname "$src")" && pwd)
        src=$(readlink "$src")
        [[ $src != /* ]] && src=$dir/$src
    done
    cd -P "$(dirname "$src")" && pwd
}

# Capture repository root and default config so other scripts can reference them.
readonly __PIPELINE_COMMON_DIR=$(_current_script_dir)
: "${PIPELINE_ROOT:=$(cd "${__PIPELINE_COMMON_DIR}/.." && pwd)}"
readonly PIPELINE_ROOT

CONFIG_FILE_DEFAULT="${PIPELINE_ROOT}/config/pipeline.env"

# Timestamped logging helpers for consistent messaging.
log() {
    local level=$1
    shift
    printf "[%s] %s %s\n" "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$level" "$*"
}

info() { log INFO "$*"; }
warn() { log WARN "$*" >&2; }
error() { log ERROR "$*" >&2; }

# Emit an error and stop execution.
die() {
    error "$*"
    exit 1
}

# Verify that required executables are available before running a stage.
require_tools() {
    local missing=()
    local tool
    for tool in "$@"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing+=("$tool")
        fi
    done
    if ((${#missing[@]})); then
        die "Missing required tools: ${missing[*]}"
    fi
}

# Convert relative paths to absolute paths rooted at the repository.
resolve_path() {
    local path=$1
    if [[ -z $path ]]; then
        printf '\n'
        return 0
    fi
    if [[ $path == /* ]]; then
        printf '%s\n' "${path%/}"
    else
        printf '%s\n' "${PIPELINE_ROOT}/${path%/}"
    fi
}

# Load pipeline configuration variables if the file exists; fall back to defaults.
load_config() {
    local config_file=${1:-$CONFIG_FILE_DEFAULT}
    if [[ -f $config_file ]]; then
        info "Loading config from ${config_file}"
        # shellcheck disable=SC1090
        source "$config_file"
    else
        warn "Config file ${config_file} not found; using defaults"
    fi
}

# Ensure a directory exists, creating it if needed.
ensure_directory() {
    local dir=$1
    [[ -d $dir ]] || mkdir -p "$dir"
}

# Resolve an executable path, allowing either absolute paths or PATH lookups.
resolve_executable() {
    local exe=$1
    if [[ -z $exe ]]; then
        die "Executable path not provided"
    fi
    if [[ $exe == */* ]]; then
        [[ -x $exe ]] || die "Executable not found or not runnable: $exe"
        printf '%s\n' "$exe"
    else
        local resolved
        resolved=$(command -v "$exe" 2>/dev/null) || die "Executable not found on PATH: $exe"
        printf '%s\n' "$resolved"
    fi
}