#!/bin/bash
# ==============================================================================
# ProbeDesign — Complete Setup Script
# ==============================================================================
#
# Sets up the ProbeDesign environment, installs all packages, builds/downloads
# bowtie indices, and runs validation tests.
#
# Supports: micromamba · mamba · conda  (auto-detected in that order)
# Platforms: macOS (Intel & Apple Silicon), Linux (x86_64 & aarch64)
# Windows:   use WSL2 to run this script (see README.md)
#
# Usage:
#   chmod +x setup_all.sh
#   ./setup_all.sh               # full setup
#   ./setup_all.sh --skip-genome # skip ~6.7 GB genome downloads
#   ./setup_all.sh --help
#
# Prerequisites:
#   - micromamba OR mamba OR conda installed
#   - ~8 GB free disk  (84 MB pseudogene indices + 6.7 GB genome indices)
# ==============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────

ENV_NAME="probedesign"

# Get repository root (directory containing this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
ENV_FILE="$REPO_ROOT/environment.yml"
INDEX_DIR="$REPO_ROOT/bowtie_indexes"
FASTA_DIR="$REPO_ROOT/probedesign/pseudogeneDBs"

# Colours
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Flags
SKIP_GENOME=false

for arg in "$@"; do
    case "$arg" in
        --skip-genome) SKIP_GENOME=true ;;
        --help|-h)
            echo "Usage: $0 [--skip-genome] [--help]"
            echo ""
            echo "  --skip-genome  Skip downloading genome indices (~6.7 GB)"
            echo "                 Use this if you only need pseudogene masking"
            echo "  --help         Show this message"
            exit 0
            ;;
    esac
done

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 0 — Detect package manager (micromamba > mamba > conda)
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo -e "${CYAN}Detecting conda package manager…${NC}"

CONDA_CMD=""
if command -v micromamba &>/dev/null; then
    CONDA_CMD="micromamba"
elif command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
else
    echo -e "${RED}ERROR: No conda package manager found.${NC}"
    echo ""
    echo "Please install one of the following:"
    echo "  • Miniforge (recommended — includes mamba + conda):"
    echo "      https://github.com/conda-forge/miniforge/releases/latest"
    echo "  • Anaconda / Miniconda:"
    echo "      https://docs.anaconda.com/free/miniconda/"
    exit 1
fi

echo -e "${GREEN}Found: ${CONDA_CMD} ($(${CONDA_CMD} --version 2>&1 | head -1))${NC}"
echo ""

# Helper: run a command inside the probedesign environment without needing
# interactive shell activation.
run_in_env() {
    "$CONDA_CMD" run -n "$ENV_NAME" "$@"
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Create environment from environment.yml
# ══════════════════════════════════════════════════════════════════════════════

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  1/5  Create '${ENV_NAME}' environment${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""

if [ ! -f "$ENV_FILE" ]; then
    echo -e "${RED}ERROR: environment.yml not found at ${ENV_FILE}${NC}"
    echo "Make sure you are running this script from the repository root:"
    echo "  ./setup_all.sh"
    exit 1
fi

# Check if environment already exists
ENV_EXISTS=false
if "$CONDA_CMD" env list 2>/dev/null | grep -qE "^${ENV_NAME}[[:space:]]"; then
    ENV_EXISTS=true
fi

if $ENV_EXISTS; then
    echo -e "${YELLOW}Environment '${ENV_NAME}' already exists.${NC}"
    echo -e "${YELLOW}Updating it from environment.yml (missing packages will be added)…${NC}"
    echo ""
    # Update existing environment
    if [ "$CONDA_CMD" = "micromamba" ]; then
        micromamba env update -n "$ENV_NAME" -f "$ENV_FILE" --prune -y
    else
        "$CONDA_CMD" env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
    fi
else
    echo "Creating environment '${ENV_NAME}' from ${ENV_FILE}…"
    echo "(This downloads ~500 MB of packages — a few minutes on first run)"
    echo ""
    if [ "$CONDA_CMD" = "micromamba" ]; then
        micromamba env create -f "$ENV_FILE" -y
    else
        "$CONDA_CMD" env create -f "$ENV_FILE"
    fi
fi

echo ""
echo -e "${GREEN}Python:  $(run_in_env python --version)${NC}"
echo -e "${GREEN}Bowtie:  $(run_in_env bowtie --version 2>&1 | head -1)${NC}"
echo -e "${GREEN}Streamlit: $(run_in_env python -c 'import streamlit; print(streamlit.__version__)')${NC}"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Install probedesign package (editable)
# ══════════════════════════════════════════════════════════════════════════════

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  2/5  Install probedesign package${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""

cd "$REPO_ROOT"
run_in_env pip install -e . --quiet

echo -e "${GREEN}probedesign $(run_in_env probedesign --version 2>&1 | tail -1)${NC}"
echo ""

# Explicit safety net for the TransQuant tab's pip dependency: some
# conda/mamba versions don't reliably re-resolve git+ URLs in the
# environment.yml "pip:" section on `env update --prune` for existing
# environments. This re-installs/updates it directly (idempotent —
# no-op if already at the pinned tag).
run_in_env pip install --quiet "git+https://github.com/JeffLeeLab/TransQuant-probe-weight.git@v0.2.0"

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Build pseudogene indices (bowtie-build from shipped FASTAs)
# ══════════════════════════════════════════════════════════════════════════════

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  3/5  Build pseudogene indices${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""

mkdir -p "$INDEX_DIR"

build_pseudogene_index() {
    local species="$1"
    local fasta="$2"
    local prefix="$3"

    # Skip if all 6 index files already present
    if ls "$INDEX_DIR/${prefix}".*.ebwt 1>/dev/null 2>&1; then
        local count
        count=$(ls "$INDEX_DIR/${prefix}".*.ebwt | wc -l | tr -d ' ')
        if [ "$count" -eq 6 ]; then
            echo -e "  ${species}: ${YELLOW}already built (skipping)${NC}"
            return 0
        fi
    fi

    if [ ! -f "$fasta" ]; then
        echo -e "  ${species}: ${RED}FASTA not found at ${fasta} — skipping${NC}"
        return 1
    fi

    local num_seqs
    num_seqs=$(grep -c "^>" "$fasta")
    echo -n "  ${species}: ${num_seqs} sequences … "

    local start_ts
    start_ts=$(date +%s)
    run_in_env bowtie-build "$fasta" "$INDEX_DIR/$prefix" >/dev/null 2>&1
    local elapsed=$(( $(date +%s) - start_ts ))

    echo -e "${GREEN}done (${elapsed}s)${NC}"
}

build_pseudogene_index "Human"       "$FASTA_DIR/human.fasta"       "humanPseudo"
build_pseudogene_index "Mouse"       "$FASTA_DIR/mouse.fasta"       "mousePseudo"
build_pseudogene_index "Drosophila"  "$FASTA_DIR/drosophila.fasta"  "drosophilaPseudo"

echo ""

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Download genome indices (pre-built Bowtie 2 .bt2 from AWS)
# ══════════════════════════════════════════════════════════════════════════════
#
# KEY FINDING: Bowtie 1 (v1.3+) reads Bowtie 2 .bt2 index files directly.
# Pre-built indices from AWS take ~12 min to download vs 1–3 hrs to build.
# Source: https://benlangmead.github.io/aws-indexes/bowtie
# ══════════════════════════════════════════════════════════════════════════════

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  4/5  Download genome indices${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""

if $SKIP_GENOME; then
    echo -e "${YELLOW}Skipped (--skip-genome). Genome masking will not be available.${NC}"
    echo ""
else

    download_genome_index() {
        local species="$1"
        local url="$2"
        local zip_name="$3"
        local source_prefix="$4"
        local target_prefix="$5"

        # Skip if all 6 index files already present
        if ls "$INDEX_DIR/${target_prefix}".*.bt2 1>/dev/null 2>&1; then
            local count
            count=$(ls "$INDEX_DIR/${target_prefix}".*.bt2 | wc -l | tr -d ' ')
            if [ "$count" -eq 6 ]; then
                echo -e "  ${species}: ${YELLOW}already downloaded (skipping)${NC}"
                return 0
            fi
        fi

        echo -n "  ${species}: downloading … "
        local start_ts
        start_ts=$(date +%s)
        curl -sS -L -o "$INDEX_DIR/$zip_name" "$url"
        local dl_elapsed=$(( $(date +%s) - start_ts ))
        local dl_size
        dl_size=$(du -h "$INDEX_DIR/$zip_name" | cut -f1)
        echo -n "${dl_size} in ${dl_elapsed}s … extracting … "

        cd "$INDEX_DIR"
        unzip -q -o "$zip_name"

        # Some AWS zips extract into a subdirectory — move files up
        if [ -d "$source_prefix" ]; then
            mv "$source_prefix"/*.bt2 . 2>/dev/null || true
            rm -rf "$source_prefix"
        fi

        # Rename files if the source prefix differs from the expected target prefix
        if [ "$source_prefix" != "$target_prefix" ]; then
            for file in "${source_prefix}".*.bt2; do
                [ -f "$file" ] && mv "$file" "${file/$source_prefix/$target_prefix}"
            done
        fi

        rm -f "$zip_name"
        cd "$REPO_ROOT"

        local elapsed=$(( $(date +%s) - start_ts ))
        echo -e "${GREEN}done (${elapsed}s total)${NC}"
    }

    # Human GRCh38  (~3.5 GB)
    download_genome_index "Human GRCh38" \
        "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip" \
        "GRCh38_noalt_as.zip" \
        "GRCh38_noalt_as" \
        "GCA_000001405.15_GRCh38_no_alt_analysis_set"

    # Mouse mm10  (~3.1 GB)
    download_genome_index "Mouse mm10" \
        "https://genome-idx.s3.amazonaws.com/bt/mm10.zip" \
        "mm10.zip" \
        "mm10" \
        "mm10"

    # Drosophila BDGP6  (~176 MB)
    download_genome_index "Drosophila BDGP6" \
        "https://genome-idx.s3.amazonaws.com/bt/BDGP6.zip" \
        "BDGP6.zip" \
        "BDGP6" \
        "drosophila"

    echo ""
fi

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Validation tests
# ══════════════════════════════════════════════════════════════════════════════

echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  5/5  Validation tests${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""

cd "$REPO_ROOT"
PASS=0
FAIL=0

check_pass() { echo -e "${GREEN}PASS${NC}"; ((PASS++)) || true; }
check_fail() { echo -e "${RED}FAIL${NC}"; ((FAIL++)) || true; }
check_partial() { local msg="$1"; echo -e "${YELLOW}${msg}${NC}"; ((FAIL++)) || true; }

# ── Test A: probedesign version ───────────────────────────────────────────────
echo -n "  [A] probedesign command works … "
if run_in_env probedesign --version 2>&1 | grep -q "0.1.0"; then
    check_pass
else
    check_fail
fi

# ── Test B: bowtie version ────────────────────────────────────────────────────
echo -n "  [B] bowtie command works … "
if run_in_env bowtie --version 2>&1 | head -1 | grep -q "bowtie-align-s version"; then
    check_pass
else
    check_fail
fi

# ── Test C: CDKN1A_32 (no bowtie indices needed) ─────────────────────────────
echo -n "  [C] probe design — CDKN1A_32 (32 probes, 100% match expected) … "
TMPDIR_C=$(mktemp -d)
run_in_env probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
    --probes 32 \
    --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa \
    -o "$TMPDIR_C/test" --quiet 2>/dev/null

ACTUAL_COUNT=$(wc -l < "$TMPDIR_C/test_oligos.txt" | tr -d ' ')
if [ "$ACTUAL_COUNT" -eq 32 ]; then
    awk -F'\t' '{print $(NF-1)}' "$TMPDIR_C/test_oligos.txt" | sort > "$TMPDIR_C/actual.txt"
    awk -F'\t' '{print $(NF-1)}' test_cases/CDKN1A_32/CDKN1A_32_genomemaskoff_oligos.txt | sort > "$TMPDIR_C/expected.txt"
    if diff -q "$TMPDIR_C/actual.txt" "$TMPDIR_C/expected.txt" >/dev/null 2>&1; then
        echo -e "${GREEN}PASS (32/32 match)${NC}"; ((PASS++)) || true
    else
        MATCH=$(comm -12 "$TMPDIR_C/actual.txt" "$TMPDIR_C/expected.txt" | wc -l | tr -d ' ')
        echo -e "${RED}FAIL (${MATCH}/32 match)${NC}"; ((FAIL++)) || true
    fi
else
    echo -e "${RED}FAIL (expected 32 probes, got ${ACTUAL_COUNT})${NC}"; ((FAIL++)) || true
fi
rm -rf "$TMPDIR_C"

# ── Test D: pseudogene index files ────────────────────────────────────────────
echo -n "  [D] pseudogene index files present … "
PSEUDO_COUNT=$(ls "$INDEX_DIR"/*Pseudo*.ebwt 2>/dev/null | wc -l | tr -d ' ')
if [ "$PSEUDO_COUNT" -ge 18 ]; then
    echo -e "${GREEN}PASS (${PSEUDO_COUNT} files)${NC}"; ((PASS++)) || true
else
    check_partial "PARTIAL (${PSEUDO_COUNT}/18 files — run Section 3 again)"
fi

if ! $SKIP_GENOME; then
    # ── Test E: genome index files ────────────────────────────────────────────
    echo -n "  [E] genome index files present … "
    GENOME_COUNT=$(ls "$INDEX_DIR"/*.bt2 2>/dev/null | wc -l | tr -d ' ')
    if [ "$GENOME_COUNT" -ge 18 ]; then
        echo -e "${GREEN}PASS (${GENOME_COUNT} files)${NC}"; ((PASS++)) || true
    else
        check_partial "PARTIAL (${GENOME_COUNT}/18 files — run Section 4 again)"
    fi

    # ── Test F: bowtie queries against all indices ────────────────────────────
    echo -n "  [F] bowtie queries against all 6 indices … "
    TEST_FA=$(mktemp /tmp/probedesign_test_XXXX.fa)
    echo -e ">test\nATCGATCGATCGATCG" > "$TEST_FA"
    IDX_PASS=0
    IDX_TOTAL=6

    for idx in humanPseudo mousePseudo drosophilaPseudo \
               GCA_000001405.15_GRCh38_no_alt_analysis_set mm10 drosophila; do
        if run_in_env bowtie -f -v 0 -k 1 --quiet \
                "$INDEX_DIR/$idx" "$TEST_FA" >/dev/null 2>&1; then
            ((IDX_PASS++)) || true
        fi
    done
    rm -f "$TEST_FA"

    if [ "$IDX_PASS" -eq "$IDX_TOTAL" ]; then
        echo -e "${GREEN}PASS (${IDX_PASS}/${IDX_TOTAL})${NC}"; ((PASS++)) || true
    else
        check_partial "PARTIAL (${IDX_PASS}/${IDX_TOTAL} indices queryable)"
    fi

    # ── Test G: KRT19 with full masking ───────────────────────────────────────
    echo -n "  [G] probe design — KRT19 (6 probes, pseudogene+genome mask, 100% match) … "
    TMPDIR_G=$(mktemp -d)
    if run_in_env probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa \
            --probes 32 --pseudogene-mask --genome-mask \
            --index-dir "$INDEX_DIR" \
            -o "$TMPDIR_G/test" --quiet 2>/dev/null; then

        G_COUNT=$(wc -l < "$TMPDIR_G/test_oligos.txt" | tr -d ' ')
        if [ "$G_COUNT" -eq 6 ]; then
            awk -F'\t' '{print $(NF-1)}' "$TMPDIR_G/test_oligos.txt" | sort > "$TMPDIR_G/actual.txt"
            awk -F'\t' '{print $(NF-1)}' test_cases/KRT19_withUTRs/KRT19_withUTRs_oligos.txt | sort > "$TMPDIR_G/expected.txt"
            if diff -q "$TMPDIR_G/actual.txt" "$TMPDIR_G/expected.txt" >/dev/null 2>&1; then
                echo -e "${GREEN}PASS (6/6 match)${NC}"; ((PASS++)) || true
            else
                MATCH=$(comm -12 "$TMPDIR_G/actual.txt" "$TMPDIR_G/expected.txt" | wc -l | tr -d ' ')
                echo -e "${RED}FAIL (${MATCH}/6 match)${NC}"; ((FAIL++)) || true
            fi
        else
            echo -e "${RED}FAIL (expected 6 probes, got ${G_COUNT})${NC}"; ((FAIL++)) || true
        fi
    else
        echo -e "${RED}FAIL (command error)${NC}"; ((FAIL++)) || true
    fi
    rm -rf "$TMPDIR_G"
fi

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  Setup Complete${NC}"
echo -e "${BLUE}======================================================${NC}"
echo ""
TOTAL=$((PASS + FAIL))
echo "  Tests: ${PASS}/${TOTAL} passed"
if [ "$FAIL" -gt 0 ]; then
    echo -e "  ${YELLOW}${FAIL} check(s) did not fully pass — see output above${NC}"
else
    echo -e "  ${GREEN}All checks passed!${NC}"
fi
echo ""
echo "  ┌─ Next steps ──────────────────────────────────────────────┐"
echo "  │                                                           │"
echo "  │  Activate the environment:                                │"
echo "  │    ${CONDA_CMD} activate ${ENV_NAME}                              │"
echo "  │                                                           │"
echo "  │  Launch the web app:                                      │"
echo "  │    streamlit run streamlit_app/app.py                     │"
echo "  │                                                           │"
echo "  │  Or use the command line:                                 │"
echo "  │    probedesign design input.fa --probes 32                │"
echo "  │                                                           │"
echo "  └───────────────────────────────────────────────────────────┘"
echo ""
