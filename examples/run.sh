#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# A new directory protects the supplied example results from overwriting.
output="${1:-results}"
mkdir "$output"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=2

python3 ../cap_auto.py inputs/n2_pure.mat --output "$output"
