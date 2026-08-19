#!/usr/bin/env bash
# Run all embedding methods on the given datasets.
# Output: outputs/<dataset>/<method>/embedding.tsv.gz, one log per method.
# Usage: bash run_all.sh <dataset> [<dataset> ...]
set -uo pipefail
export BR_ATAC_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$BR_ATAC_ROOT"

DATASETS=("${@:-Weinand}")
SERIAL="${BR_SERIAL:-0}"

# Mallet threads per pycisTopic job.
export BR_MALLET_CPU="${BR_MALLET_CPU:-16}"

run_method() {                                  # <dataset> <env> <cmd...>
  local d=$1 env=$2; shift 2
  local name; name=$(basename "$1" | sed 's/^run_//; s/\.\(R\|py\)$//')
  local log="outputs/${d}.${name}.log"
  echo "   -> $name (log: $log)"
  if [ "$1" != "${1%.R}" ]; then
    conda run --no-capture-output -n "$env" Rscript "$@" "$d" > "$log" 2>&1
  else
    conda run --no-capture-output -n "$env" python "$@" "$d" > "$log" 2>&1
  fi
  local rc=$?
  [ $rc -eq 0 ] && echo "   OK  $name ($d)" || echo "   !!  FAILED $name ($d) rc=$rc -- see $log"
}

for d in "${DATASETS[@]}"; do
  echo "=== dataset: $d  ($(date +%H:%M:%S)) ==="
  # prep must complete first: every method reads its sidecars. Cheap (seconds).
  conda run --no-capture-output -n BR_ATAC_py python common/manifest.py prep "$d" || {
    echo "!! prep failed for $d -- skipping"; continue; }

  if [ "$SERIAL" = "1" ]; then
    run_method "$d" BR_ATAC_R          methods/run_signac.R
    run_method "$d" BR_ATAC_py         methods/run_tfidf_svd.py
    run_method "$d" BR_ATAC_py         methods/run_snapatac2.py
    run_method "$d" BR_ATAC_py         methods/run_peakvi.py
    run_method "$d" BR_ATAC_pycistopic methods/run_pycistopic.py
  else
    run_method "$d" BR_ATAC_R          methods/run_signac.R      &   # + seurat_rlsi
    run_method "$d" BR_ATAC_py         methods/run_tfidf_svd.py  &   # TF-IDF+SVD -> X_TF-IDF_SVD
    run_method "$d" BR_ATAC_py         methods/run_snapatac2.py  &
    run_method "$d" BR_ATAC_py         methods/run_peakvi.py     &   # only GPU-bound one
    run_method "$d" BR_ATAC_pycistopic methods/run_pycistopic.py &   # slowest (LDA)
    wait
  fi
  echo "=== done: $d  ($(date +%H:%M:%S)) ==="
done
