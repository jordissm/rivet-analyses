#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Usage
# -----------------------------
usage() {
    echo "Usage: $0 <project> [output-path]"
    echo ""
    echo "Example:"
    echo "  $0 project"
    echo "  $0 project $PWD"
    exit 1
}


# -----------------------------
# Parse arguments
# -----------------------------
if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage
fi

PROJECT="$1"
OUTPUT_PATH="${2:-$PWD}"
BUILD_DIR="$OUTPUT_PATH/build/$PROJECT"

# -----------------------------
# Validate project directory
# -----------------------------
if [[ ! -d "$PROJECT" ]]; then
    echo "Error: Project directory '$PROJECT' does not exist."
    exit 1
fi

# -----------------------------
# Create and validate build directory
# -----------------------------
echo "Creating build directory"
mkdir -p "$BUILD_DIR"
if [[ ! -d "$BUILD_DIR" ]]; then
    echo "Error: Build directory '$BUILD_DIR' could not be created."
    exit 1
fi

# -----------------------------
# Build Rivet analyses
# -----------------------------
echo "Building Rivet analyses in '$PROJECT'..."
echo "Output directory: '$BUILD_DIR'"
echo ""

shopt -s nullglob

FOUND=0

for src in "$PROJECT"/*.cc; do
    FOUND=1

    base=$(basename "$src" .cc)
    build="${BUILD_DIR}/Rivet_${base}.so"

    echo "[BUILD] $src -> $build"

    rivet-build "$build" "$src"
done

if [[ $FOUND -eq 0 ]]; then
    echo "Warning: No .cc files found in '$PROJECT'"
fi

echo ""
echo "Done."