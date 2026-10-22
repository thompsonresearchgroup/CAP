#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
: "${GAUSSIAN_EXE:=g16}"
: "${FORMCHK_EXE:=formchk}"
: "${GAUSS_SCRDIR:=/tmp/g16-cap-${USER}}"
export GAUSS_SCRDIR
mkdir -p "$GAUSS_SCRDIR"
for stem in one_center_spdfghi_pure one_center_spdfghi_cartesian \
            h2_pure h2_cartesian n2_pure n2_cartesian \
            n2_pure_translated_rotated h2_rotated he2_spdfghi_pure he2_spdfghi_cartesian h_uhf n2_spd_bruteforce; do
  "$GAUSSIAN_EXE" < "$stem.gjf" > "$stem.log"
  "$FORMCHK_EXE" "$stem.chk" "$stem.fchk"
  grep -q 'Normal termination of Gaussian 16' "$stem.log"
done
