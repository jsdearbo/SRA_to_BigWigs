#!/usr/bin/env bash
# Smoke tests for SRA_to_BigWigs pipeline — no bioinformatics tools required.
# Run from the repository root: bash tests/test_pipeline.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PASS=0
FAIL=0

pass() { echo "PASS: $1"; (( PASS++ )) || true; }
fail() { echo "FAIL: $1"; (( FAIL++ )) || true; }

# --- Help flags exit cleanly ---

"${REPO_ROOT}/launch.sh" --help > /dev/null 2>&1 \
    && pass "launch.sh --help exits 0" \
    || fail "launch.sh --help exits non-zero"

"${REPO_ROOT}/processing_pipeline.sh" --help > /dev/null 2>&1 \
    && pass "processing_pipeline.sh --help exits 0" \
    || fail "processing_pipeline.sh --help exits non-zero"

"${REPO_ROOT}/download_data.sh" --help > /dev/null 2>&1 \
    && pass "download_data.sh --help exits 0" \
    || fail "download_data.sh --help exits non-zero"

# --- Config example loads without error ---

(
    # Source in a subshell so exports don't leak
    set -euo pipefail
    # shellcheck source=config/pipeline.env.example
    source "${REPO_ROOT}/config/pipeline.env.example"
    [[ -n "${METADATA_FILE}" ]] || exit 1
    [[ -n "${BASE_DIR}" ]]      || exit 1
    [[ -n "${ALIGN_THREADS}" ]] || exit 1
) && pass "pipeline.env.example sources cleanly" \
  || fail "pipeline.env.example failed to source"

# --- Metadata CSV has expected columns ---

python3 - "${REPO_ROOT}/metadata_ex.csv" <<'EOF'
import csv, sys
required = {"CellType", "RunID", "Description"}
with open(sys.argv[1], newline="") as f:
    reader = csv.DictReader(f)
    missing = required - set(reader.fieldnames or [])
    if missing:
        print(f"Missing columns: {missing}", file=sys.stderr)
        sys.exit(1)
    rows = list(reader)
    if len(rows) < 1:
        print("No data rows found", file=sys.stderr)
        sys.exit(1)
    print(f"  Found {len(rows)} rows with required columns")
EOF
pass "metadata_ex.csv has required columns" || fail "metadata_ex.csv column check failed"

# --- metadata.csv has expected columns ---

python3 - "${REPO_ROOT}/metadata.csv" <<'EOF'
import csv, sys
required = {"CellType", "RunID", "Description"}
with open(sys.argv[1], newline="") as f:
    reader = csv.DictReader(f)
    missing = required - set(reader.fieldnames or [])
    if missing:
        print(f"Missing columns: {missing}", file=sys.stderr)
        sys.exit(1)
    rows = list(reader)
    if len(rows) < 1:
        print("No data rows found", file=sys.stderr)
        sys.exit(1)
    print(f"  Found {len(rows)} rows with required columns")
EOF
pass "metadata.csv has required columns" || fail "metadata.csv column check failed"

# --- Summary ---

echo ""
echo "Results: ${PASS} passed, ${FAIL} failed"
(( FAIL == 0 ))
