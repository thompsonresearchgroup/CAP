#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Exercise runtime-checked code, then leave the production optimized binary.
make strict tests/read_cap
python3 -m unittest discover -s tests -p 'test_*.py' -v
make -B all
