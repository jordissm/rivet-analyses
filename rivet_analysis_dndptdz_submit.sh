#!/usr/bin/env bash
set -euo pipefail

# rivet_analysis_dndptdz_submit.sh
#
# Run the Rivet analysis EHIJING_SMASH_DNDPTDZ once per SMASH "profile_p_XXXXXX",
# WITHOUT merging across profiles, and write one YODA per profile.
#
# Tailored for Rivet v3.1.7 (no --analysis-opts), so we pass configuration via env vars:
#   RIVET_METAFILE, RIVET_PT_MIN, RIVET_PT_MAX, RIVET_PT_NBINS
#
# Usage:
#   ./run_rivet_dndptdz_per_profile.sh \
#     /path/to/output/runs/smash/events/ \
#     /path/to/output/runs/ehijing/DISKinematics.meta.jsonl \
#     /path/to/Rivet_EHIJING_SMASH_2026_DNDPTDZ.so \
#     [PT_MIN] [PT_MAX] [PT_NBINS]
#
# Examples:
#   bash run_rivet_dndptdz_per_profile.sh output/runs/smash/events/ output/runs/ehijing/events/DISKinematics.meta.jsonl rivet/Rivet_EHIJING_SMASH_2026_DNDPTDZ.so 0.0 1.1 20

SMASH_DIR="${1:?Need SMASH dir (e.g. output/runs/smash/events/)}"
META_JSON="${2:?Need meta json path (e.g. output/runs/ehijing/DISKinematics.meta.jsonl)}"
RIVET_SO="${3:?Need Rivet analysis .so path}"
PT_MIN="${4:-0.0}"
PT_MAX="${5:-1.1}"
PT_NBINS="${6:-20}"

ANA="EHIJING_SMASH_DNDPTDZ"

# Resolve rivet binary (prefer PATH, else fall back to common install prefix)
RIVET_BIN="$(command -v rivet || true)"
if [[ -z "${RIVET_BIN}" ]]; then
  if [[ -x "/opt/electra/Rivet_install/bin/rivet" ]]; then
    RIVET_BIN="/opt/electra/Rivet_install/bin/rivet"
  else
    echo "ERROR: rivet not found in PATH and /opt/electra/Rivet_install/bin/rivet not executable" >&2
    exit 1
  fi
fi

if [[ ! -d "$SMASH_DIR" ]]; then
  echo "ERROR: SMASH dir not found: $SMASH_DIR" >&2
  exit 1
fi
if [[ ! -f "$META_JSON" ]]; then
  echo "ERROR: meta json not found: $META_JSON" >&2
  exit 1
fi
if [[ ! -f "$RIVET_SO" ]]; then
  echo "ERROR: Rivet .so not found: $RIVET_SO" >&2
  exit 1
fi

# Numeric validators
is_number() {
  [[ "${1:-}" =~ ^[-+]?[0-9]*([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]
}
is_int() {
  [[ "${1:-}" =~ ^[0-9]+$ ]]
}

USE_PT_AXIS=0
if [[ -n "$PT_MIN" || -n "$PT_MAX" || -n "$PT_NBINS" ]]; then
  # Must provide all three
  if [[ -z "$PT_MIN" || -z "$PT_MAX" || -z "$PT_NBINS" ]]; then
    echo "ERROR: Provide PT_MIN PT_MAX PT_NBINS together, or provide none to use defaults." >&2
    exit 1
  fi
  if ! is_number "$PT_MIN" || ! is_number "$PT_MAX"; then
    echo "ERROR: PT_MIN/PT_MAX must be numbers (got '$PT_MIN' '$PT_MAX')" >&2
    exit 1
  fi
  if ! is_int "$PT_NBINS"; then
    echo "ERROR: PT_NBINS must be a positive integer (got '$PT_NBINS')" >&2
    exit 1
  fi
  # Compare as floats using awk (portable): require PT_MAX > PT_MIN
  if awk "BEGIN{exit !($PT_MAX <= $PT_MIN)}"; then
    echo "ERROR: PT_MAX <= PT_MIN ($PT_MAX <= $PT_MIN)" >&2
    exit 1
  fi
  if [[ "$PT_NBINS" -le 0 ]]; then
    echo "ERROR: PT_NBINS must be > 0 (got '$PT_NBINS')" >&2
    exit 1
  fi
  USE_PT_AXIS=1
fi

# Output directory
OUTDIR="/workspace/output/runs/rivet/rivet_out_${ANA}"
if [[ $USE_PT_AXIS -eq 1 ]]; then
  PT_TAG="pt_${PT_MIN}_to_${PT_MAX}_nb${PT_NBINS}"
  PT_TAG="${PT_TAG//./p}"
  PT_TAG="${PT_TAG//-/m}"
  OUTDIR="${OUTDIR}_${PT_TAG}"
fi
mkdir -p "$OUTDIR"

# Find all HepMC files
mapfile -t FILES < <(find "$SMASH_DIR" -type f -name "SMASH_HepMC_particles.asciiv3" | sort)
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "ERROR: No SMASH_HepMC_particles.asciiv3 found under $SMASH_DIR" >&2
  exit 1
fi

# Extract distinct profile names (profile_p_XXXXXX)
PROFILES=$(printf "%s\n" "${FILES[@]}" \
  | sed -n 's|^.*/\(profile_p_[0-9]\{6\}\)/SMASH_HepMC_particles\.asciiv3$|\1|p' \
  | sort -u)

echo "Rivet bin: ${RIVET_BIN}"
echo "Analysis:  ${ANA}"
echo "Meta file: ${META_JSON}"
echo "Plugin:    ${RIVET_SO}"
if [[ $USE_PT_AXIS -eq 1 ]]; then
  echo "pT axis:   ${PT_NBINS} bins in [${PT_MIN}, ${PT_MAX}] GeV"
else
  echo "pT axis:   (defaults from analysis: PT_MIN=0, PT_MAX=2, PT_NBINS=40)"
fi
echo "Outdir:    ${OUTDIR}"
echo
echo "Found profiles:"
echo "$PROFILES" | sed 's|^|  - |'

ANA_PATH="$(dirname "$RIVET_SO")"

# Run once per profile so we do NOT merge across profiles
while IFS= read -r PROFILE_NAME; do
  [[ -z "$PROFILE_NAME" ]] && continue

  YODA_OUT="${OUTDIR}/${PROFILE_NAME}.yoda"

  # Collect the matching profile file across all evt_XXXXXXXX (same profile id)
  mapfile -t PF_FILES < <(find "$SMASH_DIR" -type f -path "*/${PROFILE_NAME}/SMASH_HepMC_particles.asciiv3" | sort)

  if [[ ${#PF_FILES[@]} -eq 0 ]]; then
    echo "WARNING: No files found for ${PROFILE_NAME}, skipping." >&2
    continue
  fi

  echo
  echo "=== Running Rivet for ${PROFILE_NAME} ==="
  echo "Nfiles:  ${#PF_FILES[@]}"
  echo "Output:  ${YODA_OUT}"
  echo

  # Build env for this run
  export RIVET_METAFILE="$META_JSON"

  if [[ $USE_PT_AXIS -eq 1 ]]; then
    export RIVET_PT_MIN="$PT_MIN"
    export RIVET_PT_MAX="$PT_MAX"
    export RIVET_PT_NBINS="$PT_NBINS"
  else
    unset RIVET_PT_MIN || true
    unset RIVET_PT_MAX || true
    unset RIVET_PT_NBINS || true
  fi

  # Print the exact command
  echo "RIVET_METAFILE=\"${RIVET_METAFILE}\" \\"
  if [[ $USE_PT_AXIS -eq 1 ]]; then
    echo "RIVET_PT_MIN=\"${RIVET_PT_MIN}\" \\"
    echo "RIVET_PT_MAX=\"${RIVET_PT_MAX}\" \\"
    echo "RIVET_PT_NBINS=\"${RIVET_PT_NBINS}\" \\"
  fi
  echo "\"${RIVET_BIN}\" --ignore-beams -a \"${ANA}\" --analysis-path \"${ANA_PATH}\" -o \"${YODA_OUT}\" \\"
  echo "  <${#PF_FILES[@]} input files>"
  echo

  "${RIVET_BIN}" \
    -a "${ANA}" \
    --analysis-path "${ANA_PATH}" \
    --ignore-beams \
    -o "${YODA_OUT}" \
    "${PF_FILES[@]}"

done <<< "$PROFILES"

echo
echo "Done. Outputs in: ${OUTDIR}"