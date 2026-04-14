#!/usr/bin/env bash
#
# Compare sracha output against fasterq-dump across split modes and formats.
#
# Usage: bash validation/compare_with_fasterq_dump.sh
#
# Requires: SRR28588231.sra in crates/sracha-core/tests/fixtures/
#           (auto-downloaded by the Rust integration tests)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SRA_FILE="$ROOT_DIR/crates/sracha-core/tests/fixtures/SRR28588231.sra"
ACCESSION="SRR28588231"
SRACHA="$ROOT_DIR/target/release/sracha"
SRATOOLS_DIR="$SCRIPT_DIR/sra-tools"
COMPARE_PY="$SCRIPT_DIR/compare_fastq.py"

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
declare -a RESULTS=()

# ---------- colors ----------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BOLD='\033[1m'
RESET='\033[0m'

# ---------- helpers ----------

log()  { echo -e "${BOLD}==> $1${RESET}"; }
pass() { echo -e "  ${GREEN}PASS${RESET} $1"; PASS_COUNT=$((PASS_COUNT + 1)); RESULTS+=("PASS: $1"); }
fail() { echo -e "  ${RED}FAIL${RESET} $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); RESULTS+=("FAIL: $1"); }
skip() { echo -e "  ${YELLOW}SKIP${RESET} $1"; SKIP_COUNT=$((SKIP_COUNT + 1)); RESULTS+=("SKIP: $1"); }

# Compare two files by md5. If they differ, run compare_fastq.py for details.
# If md5 differs but seq/qual are identical (defline-only diff), still counts as PASS.
compare_files() {
    local file_a="$1"
    local file_b="$2"
    local label="$3"

    if [[ ! -f "$file_a" ]]; then
        fail "$label — sracha output missing: $(basename "$file_a")"
        return
    fi
    if [[ ! -f "$file_b" ]]; then
        fail "$label — fasterq-dump output missing: $(basename "$file_b")"
        return
    fi

    local md5_a md5_b
    md5_a=$(md5sum "$file_a" | awk '{print $1}')
    md5_b=$(md5sum "$file_b" | awk '{print $1}')

    if [[ "$md5_a" == "$md5_b" ]]; then
        local size
        size=$(stat --printf='%s' "$file_a")
        pass "$label — byte-identical (md5=$md5_a, size=$size)"
        return
    fi

    echo -e "  ${YELLOW}md5 mismatch${RESET}: sracha=$md5_a  fasterq-dump=$md5_b"

    # For FASTQ files, check if seq/qual match (defline-only diff is OK)
    local ext="${file_a##*.}"
    if [[ "$ext" == "fastq" ]]; then
        echo "  Checking sequence/quality identity with compare_fastq.py..."
        if python3 "$COMPARE_PY" "$file_a" "$file_b"; then
            pass "$label — content-identical (deflines differ)"
            echo "  Showing defline difference (first record):"
            echo "    sracha:      $(head -1 "$file_a")"
            echo "    fasterq-dump: $(head -1 "$file_b")"
        else
            fail "$label — seq/qual mismatch"
        fi
    elif [[ "$ext" == "fasta" ]]; then
        # For FASTA: check if diffs are only N-masking (fasterq-dump replaces
        # low-quality bases with N in FASTA mode, sracha preserves originals)
        local n_mask_result
        n_mask_result=$(paste <(grep -v '^>' "$file_a") <(grep -v '^>' "$file_b") | python3 -c "
import sys
total = n_mask = other = 0
n_positions = 0
for line in sys.stdin:
    parts = line.rstrip().split('\t')
    if len(parts) != 2: continue
    a, b = parts
    total += 1
    if a == b: continue
    all_n = True
    for ca, cb in zip(a, b):
        if ca != cb:
            if cb == 'N': n_positions += 1
            else: all_n = False
    if all_n: n_mask += 1
    else: other += 1
print(f'{total} {n_mask} {other} {n_positions}')
")
        local total n_mask other n_positions
        read -r total n_mask other n_positions <<< "$n_mask_result"

        if [[ "$other" -eq 0 ]]; then
            if [[ "$n_mask" -eq 0 ]]; then
                # Only header diffs
                pass "$label — sequences identical (headers differ)"
            else
                pass "$label — N-masking only ($n_mask/$total seqs, $n_positions positions; fasterq-dump masks Q2 bases)"
            fi
            # Show header diff if present
            local hdr_a hdr_b
            hdr_a=$(head -1 "$file_a")
            hdr_b=$(head -1 "$file_b")
            if [[ "$hdr_a" != "$hdr_b" ]]; then
                echo "    header diff: sracha='$hdr_a' vs fasterq='$hdr_b'"
            fi
        else
            echo "  FASTA has $other non-N-masking diffs — showing first diff:"
            diff <(head -20 "$file_a") <(head -20 "$file_b") || true
            fail "$label — sequence mismatch ($other diffs beyond N-masking)"
        fi
    else
        fail "$label — unknown extension .$ext"
    fi
}

# ---------- preflight ----------

log "Checking prerequisites..."

if [[ ! -f "$SRA_FILE" ]]; then
    echo "ERROR: Test fixture not found: $SRA_FILE"
    echo "Run the integration tests first: cargo test -p sracha-core -- --ignored"
    exit 1
fi

if [[ ! -x "$SRACHA" ]]; then
    echo "ERROR: sracha binary not found at $SRACHA"
    echo "Build with: cargo build --release -p sracha"
    exit 1
fi

echo "  sracha: $($SRACHA --version)"
echo "  SRA file: $SRA_FILE ($(du -h "$SRA_FILE" | awk '{print $1}'))"

# ---------- install sra-tools if needed ----------

FASTERQ_DUMP=""
# check if already installed locally
if compgen -G "$SRATOOLS_DIR/sratoolkit.*/bin/fasterq-dump" > /dev/null 2>&1; then
    FASTERQ_DUMP=$(ls "$SRATOOLS_DIR"/sratoolkit.*/bin/fasterq-dump | head -1)
elif command -v fasterq-dump &>/dev/null; then
    FASTERQ_DUMP=$(command -v fasterq-dump)
fi

if [[ -z "$FASTERQ_DUMP" ]]; then
    log "Installing sra-tools to $SRATOOLS_DIR..."
    mkdir -p "$SRATOOLS_DIR"

    TARBALL_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
    TARBALL="$SRATOOLS_DIR/sratoolkit.tar.gz"

    echo "  Downloading from $TARBALL_URL ..."
    curl -fSL -o "$TARBALL" "$TARBALL_URL"
    echo "  Extracting..."
    tar -xzf "$TARBALL" -C "$SRATOOLS_DIR"
    rm -f "$TARBALL"

    FASTERQ_DUMP=$(ls "$SRATOOLS_DIR"/sratoolkit.*/bin/fasterq-dump | head -1)
    if [[ ! -x "$FASTERQ_DUMP" ]]; then
        echo "ERROR: fasterq-dump not found after extraction"
        ls -la "$SRATOOLS_DIR"/
        exit 1
    fi
fi

echo "  fasterq-dump: $("$FASTERQ_DUMP" --version 2>&1 | head -1 || echo 'unknown version')"

# ---------- setup temp dir ----------

TMPDIR_BASE=$(mktemp -d "${TMPDIR:-/tmp}/sracha-compare.XXXXXX")
trap 'rm -rf "$TMPDIR_BASE"' EXIT
echo "  Temp dir: $TMPDIR_BASE"
echo

# =====================================================================
# TEST 1: split-3 (paired-end, default)
# =====================================================================
log "Test 1: split-3 (FASTQ, paired-end)"

SRACHA_OUT="$TMPDIR_BASE/test1_sracha"
FASTERQ_OUT="$TMPDIR_BASE/test1_fasterq"
mkdir -p "$SRACHA_OUT" "$FASTERQ_OUT"

"$SRACHA" fastq "$SRA_FILE" --split split-3 --no-gzip -O "$SRACHA_OUT" -f --no-progress 2>&1 | tail -2
"$FASTERQ_DUMP" "$SRA_FILE" --split-3 -O "$FASTERQ_OUT" -f 2>&1 | tail -2

echo "  sracha files:      $(ls "$SRACHA_OUT"/)"
echo "  fasterq-dump files: $(ls "$FASTERQ_OUT"/)"

compare_files "$SRACHA_OUT/${ACCESSION}_1.fastq" "$FASTERQ_OUT/${ACCESSION}_1.fastq" "split-3 read 1"
compare_files "$SRACHA_OUT/${ACCESSION}_2.fastq" "$FASTERQ_OUT/${ACCESSION}_2.fastq" "split-3 read 2"
echo

# =====================================================================
# TEST 2: split-spot
# =====================================================================
log "Test 2: split-spot (FASTQ)"

SRACHA_OUT="$TMPDIR_BASE/test2_sracha"
FASTERQ_OUT="$TMPDIR_BASE/test2_fasterq"
mkdir -p "$SRACHA_OUT" "$FASTERQ_OUT"

"$SRACHA" fastq "$SRA_FILE" --split split-spot --no-gzip -O "$SRACHA_OUT" -f --no-progress 2>&1 | tail -2
"$FASTERQ_DUMP" "$SRA_FILE" --split-spot -O "$FASTERQ_OUT" -f 2>&1 | tail -2

echo "  sracha files:      $(ls "$SRACHA_OUT"/)"
echo "  fasterq-dump files: $(ls "$FASTERQ_OUT"/)"

compare_files "$SRACHA_OUT/${ACCESSION}.fastq" "$FASTERQ_OUT/${ACCESSION}.fastq" "split-spot"
echo

# =====================================================================
# TEST 3: split-files
# =====================================================================
log "Test 3: split-files (FASTQ)"

SRACHA_OUT="$TMPDIR_BASE/test3_sracha"
FASTERQ_OUT="$TMPDIR_BASE/test3_fasterq"
mkdir -p "$SRACHA_OUT" "$FASTERQ_OUT"

"$SRACHA" fastq "$SRA_FILE" --split split-files --no-gzip -O "$SRACHA_OUT" -f --no-progress 2>&1 | tail -2
"$FASTERQ_DUMP" "$SRA_FILE" --split-files -O "$FASTERQ_OUT" -f 2>&1 | tail -2

echo "  sracha files:      $(ls "$SRACHA_OUT"/)"
echo "  fasterq-dump files: $(ls "$FASTERQ_OUT"/)"

# split-files for 2-read paired data should produce _1 and _2
for suffix in _1 _2; do
    sf="$SRACHA_OUT/${ACCESSION}${suffix}.fastq"
    ff="$FASTERQ_OUT/${ACCESSION}${suffix}.fastq"
    if [[ -f "$sf" && -f "$ff" ]]; then
        compare_files "$sf" "$ff" "split-files ${suffix}"
    elif [[ -f "$sf" ]]; then
        fail "split-files ${suffix} — fasterq-dump did not produce ${suffix}.fastq"
    elif [[ -f "$ff" ]]; then
        fail "split-files ${suffix} — sracha did not produce ${suffix}.fastq"
    fi
done
echo

# =====================================================================
# TEST 4: FASTA split-3
# =====================================================================
log "Test 4: split-3 (FASTA, paired-end)"

SRACHA_OUT="$TMPDIR_BASE/test4_sracha"
FASTERQ_OUT="$TMPDIR_BASE/test4_fasterq"
mkdir -p "$SRACHA_OUT" "$FASTERQ_OUT"

"$SRACHA" fastq "$SRA_FILE" --split split-3 --no-gzip --fasta -O "$SRACHA_OUT" -f --no-progress 2>&1 | tail -2
"$FASTERQ_DUMP" "$SRA_FILE" --split-3 --fasta -O "$FASTERQ_OUT" -f 2>&1 | tail -2

echo "  sracha files:      $(ls "$SRACHA_OUT"/)"
echo "  fasterq-dump files: $(ls "$FASTERQ_OUT"/)"

compare_files "$SRACHA_OUT/${ACCESSION}_1.fasta" "$FASTERQ_OUT/${ACCESSION}_1.fasta" "FASTA split-3 read 1"
compare_files "$SRACHA_OUT/${ACCESSION}_2.fasta" "$FASTERQ_OUT/${ACCESSION}_2.fasta" "FASTA split-3 read 2"
echo

# =====================================================================
# TEST 5: FASTA split-spot
# =====================================================================
log "Test 5: split-spot (FASTA)"

SRACHA_OUT="$TMPDIR_BASE/test5_sracha"
FASTERQ_OUT="$TMPDIR_BASE/test5_fasterq"
mkdir -p "$SRACHA_OUT" "$FASTERQ_OUT"

"$SRACHA" fastq "$SRA_FILE" --split split-spot --no-gzip --fasta -O "$SRACHA_OUT" -f --no-progress 2>&1 | tail -2
"$FASTERQ_DUMP" "$SRA_FILE" --split-spot --fasta -O "$FASTERQ_OUT" -f 2>&1 | tail -2

echo "  sracha files:      $(ls "$SRACHA_OUT"/)"
echo "  fasterq-dump files: $(ls "$FASTERQ_OUT"/)"

compare_files "$SRACHA_OUT/${ACCESSION}.fasta" "$FASTERQ_OUT/${ACCESSION}.fasta" "FASTA split-spot"
echo

# =====================================================================
# SUMMARY
# =====================================================================
echo
echo -e "${BOLD}======================================${RESET}"
echo -e "${BOLD}  SUMMARY${RESET}"
echo -e "${BOLD}======================================${RESET}"
for r in "${RESULTS[@]}"; do
    if [[ "$r" == PASS* ]]; then
        echo -e "  ${GREEN}$r${RESET}"
    elif [[ "$r" == FAIL* ]]; then
        echo -e "  ${RED}$r${RESET}"
    else
        echo -e "  ${YELLOW}$r${RESET}"
    fi
done
echo
echo -e "  Passed: ${GREEN}${PASS_COUNT}${RESET}  Failed: ${RED}${FAIL_COUNT}${RESET}  Skipped: ${YELLOW}${SKIP_COUNT}${RESET}"
echo

if [[ "$FAIL_COUNT" -gt 0 ]]; then
    echo -e "${RED}SOME TESTS FAILED${RESET}"
    exit 1
else
    echo -e "${GREEN}ALL TESTS PASSED${RESET}"
    exit 0
fi
