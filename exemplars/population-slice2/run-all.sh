#!/usr/bin/env bash
#
# run-all.sh -- the whole of WS-E slice 2, end to end.
#
#   ./fetch.sh         # 392 files, 7.21 GB, md5-verified, resumable
#   bash run-all.sh    # ~80 min on a 2023 laptop, ~4 GB peak RSS
#
# Each stage is idempotent and skips work it has already done, so an
# interrupted run can simply be repeated. Delete the corresponding file under
# data/derived/ to force a stage to recompute.
#
#   01-glm.R          per-subject per-run GLM        -> data/derived/*_betas.rds
#   02-transports.R   group grid, P^A, P^F, section  -> data/derived/transports.rds
#                     7.5 data-free diagnostics         results/...-transport-diagnostics.csv
#   03-population.R   two population fits + the      -> data/derived/population.rds
#                     identity acceptances              results/...-receipts.csv
#   04-eta-transport.R  V^C, V^W, eta and its null   -> results/...-eta.csv
#                     band, ACROSS TASK                 results/...-eta-null.csv
#   05-eta-across-run.R the same eta ACROSS RUN --   -> results/...-eta-across-run.csv
#                     the second independence axis      results/...-eta-across-run-null.csv
#
# 03, 04 and 05 must run in this order: 03 writes the receipts CSV from
# scratch and 04 and 05 merge their rows into it.
#
# Environment respected by the stages:
#   SLICE2_DIR   this directory (inferred from the script path by default)
#   SUBJECTS     comma/space separated subset, 01 only
#   ETA_DRAWS    null band draw count, 04 only (200 by default; set to 5 for a
#                smoke run that finishes in under a minute)

set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export SLICE2_DIR="$HERE"

if [[ ! -d "$HERE/data/derivatives/fmriprep" ]]; then
  printf 'run-all.sh: no data under %s/data -- run ./fetch.sh first\n' "$HERE" >&2
  exit 1
fi

for stage in 01-glm.R 02-transports.R 03-population.R 04-eta-transport.R \
             05-eta-across-run.R; do
  printf '\n===== %s =====\n' "$stage"
  Rscript "$HERE/$stage"
done

printf '\n===== done =====\n'
ls -la "$HERE/results"
