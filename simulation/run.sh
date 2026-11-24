#!/bin/zsh
# Run a Python script with single-threaded BLAS (each PyMC fit already uses 4 chains in parallel).
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec python -B -u "$@"
